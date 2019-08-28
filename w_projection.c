
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <limits.h>
#include <string.h>
#include <stdbool.h>

#include "constants.h"
#include "complex.h"
#include "fft.h"
#include "window.h"
#include "utility.h"
#include "w_projection.h"

void generate_w_projection_kernels(void)
{
    if(SINGLE_PRECISION)
        printf(">>> INFO: Generating W-Projection kernels using single precision...\n");
    else
        printf(">>> INFO: Generating W-Projection kernels using double precision...\n");
    
    int number_w_planes_to_create = 5;
    int number_w_planes = 339;
    int grid_size = 18000;
    int min_support = 4;
    int max_support = 44;
    int texture_size = 512; // MUST BE POWER OF TWO
  	
    size_t max_bytes_per_plane = 20 * 1024 * 1024; // 12MB
    
    PREC max_uvw  = 7083.386050;
    PREC w_scale = PREC_POW(number_w_planes - 1, 2.0) / max_uvw;
    PREC cell_size = 6.39708380288949E-06; //6.39954059065e-06; <= we suspect new cell size
    PREC w_to_max_support_ratio = (max_support - min_support) / max_uvw;
    PREC fov = cell_size * grid_size;

    // Ensure texture size is a power of two
    if(!is_power_of_two(texture_size))
    {
        printf(">>> ERROR: Texture size must be a power of two, exiting...\n");
        exit(EXIT_FAILURE);
    }

    // Calculate convolution kernel memory requirements
    size_t max_mem_bytes = MIN(max_bytes_per_plane * number_w_planes, get_total_ram_capacity());
    PREC max_conv_size = PREC_SQRT(max_mem_bytes / (16.0 * number_w_planes));
    printf("Max conv size: %f\n", max_conv_size);
    int nearest = get_next_pow_2((unsigned int) 2 * (int) (max_conv_size / 2.0));
    printf("Nearest: %d\n", nearest);
    int conv_size = nearest;
    printf("Conv size: %d\n", conv_size);
    int conv_half_size = conv_size / 2;
    printf("Conv half size: %d\n", conv_half_size);

    // Ensure max convolution size is larger than kernel texture size
    if(texture_size > max_conv_size)
    {
    	printf(">>> ERROR: Texture size must be less than or equal to the maximum convolution size (resolution), exiting...\n");
    	exit(EXIT_FAILURE);
    }
   
    PREC max_l = PREC_SIN(0.5 * fov);
    
    // Allocation of memory
    Complex *kernels = calloc(number_w_planes * texture_size * texture_size, sizeof(Complex));
    Complex *screen = calloc(conv_size * conv_size, sizeof(Complex));
    Complex *interpolated = calloc(texture_size * texture_size, sizeof(Complex));
    PREC* taper = calloc(conv_size, sizeof(PREC));
    
    unsigned int plane_to_save = 0;
    int start_plane = 0; // number_w_planes - number_w_planes_to_create;
    int end_plane = 5;

    printf(">>> UPDATE: Creating w projection kernels...\n");
    for(int iw = start_plane; iw < end_plane; ++iw)
    {
        printf(">>> UPDATE: Creating kernel %d\n", iw);
        memset(screen, 0, conv_size * conv_size * sizeof(Complex));
        memset(interpolated, 0, texture_size * texture_size * sizeof(Complex));
        memset(taper, 0, conv_size * sizeof(PREC));
        
        // Populate taper based on current kernel support
        PREC w = iw * iw / w_scale;
        int support = (int) PREC_ROUND(calculate_support(w, min_support, w_to_max_support_ratio));
        int full_support = support * 2 + 1;
        populate_centered_prolate(taper, conv_size, full_support);

        // Generate screen
        PREC sampling = ((2.0 * max_l) / grid_size) * ((PREC) grid_size / (PREC) full_support);
        generate_phase_screen(iw, conv_size, sampling, w_scale, taper, screen);

		if(iw == plane_to_save)        
        	save_kernel_to_file("../kernels/phase_screen.csv", screen, conv_size);

        printf(">>> UPDATE: Executing Fourier Transform...\n");
        // FFT
        fft_2d(screen, conv_size);

		if(iw == plane_to_save)        
        	save_kernel_to_file("../kernels/post-FFT.csv", screen, conv_size);

        // Interpolate kernel down to texture
        interpolate_kernel(screen, interpolated, conv_size, texture_size);

		if(iw == plane_to_save)        
        	save_kernel_to_file("../kernels/interpolated.csv", interpolated, texture_size);

        // Normalize texture kernel for scaling on GPU
        normalize_kernel(interpolated, texture_size, support);

		if(iw == plane_to_save)        
        	save_kernel_to_file("../kernels/normalized.csv", interpolated, texture_size);
        
        printf(">>> UPDATE: Storing useful quadrant for further processing...\n");
        // Clip
        for(int row_index = 0; row_index < texture_size; ++row_index)
            for(int col_index = 0; col_index < texture_size; ++col_index)
            {
                int offset = iw * texture_size * texture_size;
                int k_index = offset + row_index * texture_size + col_index;
                kernels[k_index] = interpolated[row_index * texture_size + col_index];
            }
    }
    
    free(taper);
    free(screen);
    free(interpolated);
    
    printf(">>> UPDATE: Extracting support kernel from quadrants...\n");
    FILE *kernel_real_file = fopen("../kernels/w-proj_kernels_real.csv", "w");
    FILE *kernel_imag_file = fopen("../kernels/w-proj_kernels_imag.csv", "w");
    FILE *support_file = fopen("../kernels/w-proj_supports.csv", "w");

    for(int iw = 0; iw < number_w_planes; ++iw)
    {
        PREC w = iw * iw / w_scale;
        int kernel_offset = iw * texture_size * texture_size;
        PREC support = calculate_support(w, min_support, w_to_max_support_ratio);

        fprintf(support_file, "%d\n", (int) PREC_ROUND(support));
        
        for(int row = 0; row < 1; ++row)
        {
            for(int col = 0; col < texture_size; ++col)
            {
                int plane_index = kernel_offset + row * texture_size + col;
                
                fprintf(kernel_real_file, "%.10f ", kernels[plane_index].real);
                fprintf(kernel_imag_file, "%.10f ", kernels[plane_index].imag);
            }
        }
        
        fprintf(kernel_real_file, "\n");
        fprintf(kernel_imag_file, "\n");
    }
    
    fclose(support_file);
    fclose(kernel_imag_file);
    fclose(kernel_real_file);
    free(kernels);
    
    printf(">>> UPDATE: W-Projection kernels successfully created, exiting...\n");
}

void generate_phase_screen(int iw, int conv_size, PREC sampling, PREC w_scale, PREC* taper, Complex *screen)
{
    PREC f = (2.0 * PI * iw * iw) / w_scale;
    printf("Kernel %d => F: %f\n", iw, f);
    int conv_size_half = conv_size / 2;
    
    for(int iy = -conv_size_half; iy < conv_size_half; ++iy)
    {
        PREC taper_y = taper[iy + conv_size_half];
        PREC m = sampling * (PREC) iy;
        PREC msq = m*m;
        int offset = (iy > -1 ? iy : (iy + conv_size)) * conv_size;
        
        for(int ix = -conv_size_half; ix < conv_size_half; ++ix)
        {
            PREC l = sampling * (PREC) ix;
            PREC rsq = l * l + msq;
            if (rsq < 1.0) {
                PREC taper_x = taper[ix + conv_size_half];
                PREC taper = taper_x * taper_y;
                int index = (offset + (ix > -1 ? ix : (ix + conv_size)));
                PREC phase = f * (PREC_SQRT(1.0 - rsq) - 1.0);
                screen[index] = (Complex) {
                    .real = taper * PREC_COS(phase),
                    .imag = taper * PREC_SIN(phase)
                };
            }
        }
    }
}

void crop_plane(Complex *plane, Complex *cropped_plane, int resolution, int support)
{
    int half_resolution = resolution/2;
    int half_support = support/2;
    
    for(int row_index = -half_support; row_index <= half_support; ++row_index)
    {
        int cropped_row = (row_index + half_support) * support;
        int plane_row = (row_index + half_resolution) * resolution;
        
        for(int col_index = -half_support; col_index <= half_support; ++col_index)
        {
            int cropped_index = cropped_row + half_support + col_index;
            int plane_index = plane_row + half_resolution + col_index;
            
            cropped_plane[cropped_index] = plane[plane_index];
        }
    }
}

PREC calculate_support(PREC w, int min_support, PREC w_max_support_ratio)
{
    // int support = (int) (PREC_ABS(w_max_support_ratio * w) + min_support);
    // return (support % 2 == 0) ? support + 1 : support;
    
    return PREC_ABS(w_max_support_ratio * w) + min_support;
}

void interpolate_kernel(Complex *src, Complex *dest, int src_support, int dest_support)
{
    // Determine "distance" between source samples in range [-1.0, 1.0]
    PREC src_stride = 2.0 / (src_support - 1.0);
    // Storage for neighbours, synthesized points
    Complex n[16], p[4];
    // Neighbours shift values (rs = row shift, cs = col shift)
    PREC rs[16], cs[16];
    PREC row_stride, col_stride = 0.0;
    
    for(int row_index = 0; row_index < dest_support; ++row_index)
    {
        // Determine relative shift for interpolation row [-1.0, 1.0]
        row_stride = calc_interpolation_stride((PREC) row_index, (PREC) dest_support);
        
        for(int col_index = 0; col_index < dest_support; ++col_index)
        {
            // Determine relative shift for interpolation col [-1.0, 1.0]
            col_stride = calc_interpolation_stride((PREC) col_index, (PREC) dest_support);
            
            // gather 16 neighbours
            get_bicubic_neighbours(row_stride, col_stride, n, rs, cs, src_support, src);
            // interpolate intermediate samples
            p[0] = interpolate_sample(n[0], n[1], n[2], n[3],
                cs[0], cs[1], cs[2], cs[3], src_stride, col_stride);
            p[1] = interpolate_sample(n[4], n[5], n[6], n[7],
                cs[4], cs[5], cs[6], cs[7], src_stride, col_stride);
            p[2] = interpolate_sample(n[8], n[9], n[10], n[11],
                cs[8], cs[9], cs[10], cs[11], src_stride, col_stride);
            p[3] = interpolate_sample(n[12], n[13], n[14], n[15],
                cs[12], cs[13], cs[14], cs[15], src_stride, col_stride);
            
            // interpolate final sample
            dest[row_index * dest_support + col_index] = interpolate_sample(p[0], p[1], p[2], p[3],
               rs[1], rs[5], rs[9], rs[13], src_stride, row_stride);
        }
    }
}

PREC calc_interpolation_stride(PREC index, PREC width)
{
    return -1.0 + ((2.0 * index + 1.0) / width);
}

int calc_relative_index(PREC x, PREC width)
{
    int offset = (x < 0.0) ? 1 : 2;
    return ((int) PREC_FLOOR(((x + 1.0) / 2.0) * (width - offset))) + 1;
}

void get_bicubic_neighbours(PREC row_stride, PREC col_stride, Complex *neighbours, PREC *row_strides,
	PREC *col_strides, int src_support, Complex *src)
{
    // determine where to start locating neighbours in source matrix
    int x = calc_relative_index(col_stride, (PREC) src_support);
    int y = calc_relative_index(row_stride, (PREC) src_support);
    // counter for active neighbour
    int neighbour_index = 0;
    // define neighbour boundaries
    int row_start = (row_stride < 0.0) ? y - 1 : y - 2;
    int row_end = (row_stride < 0.0) ? y + 3 : y + 2;
    int col_start = (col_stride < 0.0) ? x - 1 : x - 2;
    int col_end = (col_stride < 0.0) ? x + 3 : x + 2;
    int half_src_support = src_support / 2;
    // gather 16 neighbours
    for(int row_index = row_start; row_index < row_end; ++row_index)
    {   
        for(int col_index = col_start; col_index < col_end; ++col_index)
        {
            // set row and col shifts for neighbour
            row_strides[neighbour_index] = (row_stride < 0.0) ? 
            	calculate_window_stride(row_index - 1, src_support - 1) : calculate_window_stride(row_index, src_support - 1);
            col_strides[neighbour_index] = (col_stride < 0.0) ? 
            	calculate_window_stride(col_index - 1, src_support - 1) : calculate_window_stride(col_index, src_support - 1);            
            // neighbour falls out of bounds
            if(row_index < 0 || col_index < 0 || row_index >= half_src_support || col_index >= half_src_support)
                neighbours[neighbour_index] = (Complex) {.real = 0.0, .imag = 0.0};
            // neighbour exists
            else
                neighbours[neighbour_index] = src[row_index * src_support + col_index];

            ++neighbour_index;
        }
    }   
}

Complex interpolate_sample(Complex z0, Complex z1, Complex z2, Complex z3,
	PREC x0, PREC x1, PREC x2, PREC x3, PREC h, PREC x)
{
    PREC h_cube = PREC_POW(h, 3.0);
    PREC scale0 = -(x - x1) * (x - x2) * (x - x3) / (6.0 * h_cube);
    PREC scale1 =  (x - x0) * (x - x2) * (x - x3) / (2.0 * h_cube);
    PREC scale2 = -(x - x0) * (x - x1) * (x - x3) / (2.0 * h_cube);
    PREC scale3 =  (x - x0) * (x - x1) * (x - x2) / (6.0 * h_cube);
   
    Complex z = complex_scale(z0, scale0);
    z = complex_add(z, complex_scale(z1, scale1));
    z = complex_add(z, complex_scale(z2, scale2));
    z = complex_add(z, complex_scale(z3, scale3));
    
    return z;
}

void normalize_kernel(Complex *kernel, int texture_size, int kernel_support)
{    
    PREC real_sum = 0.0;

    for(int row_index = 0; row_index < texture_size; row_index++)
        for(int col_index = 0; col_index < texture_size; ++col_index)
            real_sum += kernel[row_index * texture_size + col_index].real;
    
    PREC norm_scalar = PREC_POW((PREC) texture_size / (PREC) kernel_support, 2.0) / real_sum;

    for(int row_index = 0; row_index < texture_size; ++row_index)
    {
        for(int col_index = 0; col_index < texture_size; ++col_index)
        {
            int index = row_index * texture_size + col_index;
            kernel[index] = complex_scale(kernel[index], norm_scalar);
        }
    }
}

void save_kernel_to_file(char *filename, Complex *kernel, int dimension)
{
    FILE *file = fopen(filename, "w");

    for(int r = 0; r < dimension; r++)
    {
        for(int c = 0; c < dimension; c++)
        {
            //if(r == support/2)
                fprintf(file, "%.10f ", kernel[r * dimension + c].real);
        }
        fprintf(file, "\n");
    }
    
    fclose(file);
}