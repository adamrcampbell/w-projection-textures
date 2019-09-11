
#ifndef W_PROJECTION_H
#define W_PROJECTION_H

#ifdef __cplusplus
extern "C" {
#endif
   
    void generate_w_projection_kernels(void);

    void generate_phase_screen(PREC f, PREC fov, int resolution_size, PREC* taper, Complex *screen);
    
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

	void normalize_kernel(Complex *kernel, int texture_size, int kernel_support, PREC true_support);

	void save_kernel_to_file(char *filename, Complex *kernel, int dimension);

	PREC calcSphrShift(PREC index, PREC width);

	void save_kernels_to_file(const char *real_file, const char *imag_file, 
	const Complex *kernels, const int width, const int num_kernels);
    
#ifdef __cplusplus
}
#endif

#endif /* W_PROJECTION_H */

