
#ifndef WINDOW_H
#define WINDOW_H

#ifdef __cplusplus
extern "C" {
#endif

    void populate_ps_window(PREC *window, int size);
    void populate_centered_prolate(PREC *window, int window_width, int prolate_width);
    PREC prolate_spheroidal(PREC nu);
    PREC calculate_window_stride(int index, int size);

#ifdef __cplusplus
}
#endif

#endif /* WINDOW_H */

