
#ifndef W_PROJECTION_H
#define W_PROJECTION_H

#ifdef __cplusplus
extern "C" {
#endif
   
    void generate_w_projection_kernels(void);
    void generate_phase_screen(int iw, int conv_size, PREC sampling, PREC w_scale, PREC* taper, Complex *screen);
    void normalize_kernels_by_maximum(Complex *kernels, PREC *maximums, int number_w_planes, int conv_half_size);
    void normalize_kernels_sum_of_one(Complex *kernels, int number_w_planes, int conv_half_size, int oversample);
    
    void normalize_plane(Complex *plane, int resolution);
    void crop_plane(Complex *plane, Complex *cropped_plane, int resolution, int support);
    PREC calculate_support(PREC w, int min_support, PREC w_max_support_ratio);

    void interpolate_kernel(Complex *src, Complex *dest, int src_support, int dest_support);

    PREC calc_interpolation_stride(PREC index, PREC width);

    int calc_relative_index(PREC x, PREC width);

    void get_bicubic_neighbours(PREC row_stride, PREC col_stride, Complex *neighbours, PREC *row_strides,
		PREC *col_strides, int src_support, Complex *src);

    Complex interpolate_sample(Complex z0, Complex z1, Complex z2, Complex z3,
	PREC x0, PREC x1, PREC x2, PREC x3, PREC h, PREC x);

	void normalize_kernel(Complex *kernel, int texture_size, int kernel_support);

	void save_kernel_to_file(char *filename, Complex *kernel, int dimension);
    
#ifdef __cplusplus
}
#endif

#endif /* W_PROJECTION_H */

