// $Id$
#ifndef __INC_PS_PLOT_H__
#define __INC_PS_PLOT_H__

int ps_plot(const char *string, const char *structure, const char *ssfile);

template < class BPTable >
int ps_color_plot(const char *string, const char *structure, const BPTable& bp, const char* ssfile);

#endif  //  __INC_PS_PLOT_H__
