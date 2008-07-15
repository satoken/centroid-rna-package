/* $Id$ */

#ifndef __INC_CAPI_H__
#define __INC_CAPI_H__

#ifdef __cplusplus
extern "C" {
#endif

#define CENTROID_SUCCESS 0
#define CENTROID_BAD_ALLOC 1
#define CENTROID_LENGTH_ERR 2
#define CENTROID_UNKNOWN_ERR 3

int centroidfold(const char *seq, char *str, float g);

#ifdef __cplusplus
}
#endif
#endif /*  __INC_CAPI_H__ */
