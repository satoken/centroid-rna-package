#include <iostream>
#include <vector>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cerrno>
#include <cmath>
#include <ctime>
#include "ps_plot.h"

#define  MAX(A,B)    (A)>(B)?(A):(B)
#define VERSION "1.8.3"

static short *make_pair_table(const char *structure);
static char *time_stamp(void);
static void nrerror(const char message[]);
extern "C" { int   naview_xy_coordinates(short *pair_table, float *X, float *Y); }

static char rcsid[] = "$Id$";

// from Vienna RNA package 1.8.3
static const char *RNAss_head =
"%%BeginProlog\n"
"/RNAplot 100 dict def\n"
"RNAplot begin\n"
"/fsize  14 def\n"
"/outlinecolor {0.2 setgray} bind def\n"
"/paircolor    {0.2 setgray} bind def\n"
"/seqcolor     {0   setgray} bind def\n"
"/cshow  { dup stringwidth pop -2 div fsize -3 div rmoveto show} bind def\n"
"/min { 2 copy gt { exch } if pop } bind def\n"
"/max { 2 copy lt { exch } if pop } bind def\n"
"/drawoutline {\n"
"  gsave outlinecolor newpath\n"
"  coor 0 get aload pop 0.8 0 360 arc % draw 5' circle of 1st sequence\n"
"  currentdict /cutpoint known        % check if cutpoint is defined\n"
"  {coor 0 cutpoint getinterval\n"
"   {aload pop lineto} forall         % draw outline of 1st sequence\n"
"   coor cutpoint get aload pop\n"
"   2 copy moveto 0.8 0 360 arc       % draw 5' circle of 2nd sequence\n"
"   coor cutpoint coor length cutpoint sub getinterval\n"
"   {aload pop lineto} forall}        % draw outline of 2nd sequence\n"
"  {coor {aload pop lineto} forall}   % draw outline as a whole\n"
"  ifelse\n"
"  stroke grestore\n"
"} bind def\n"
"/drawpairs {\n"
"  paircolor\n"
"  0.7 setlinewidth\n"
"  [9 3.01] 9 setdash\n"
"  newpath\n"
"  pairs {aload pop\n"
"     coor exch 1 sub get aload pop moveto\n"
"     coor exch 1 sub get aload pop lineto\n"
"  } forall\n"
"  stroke\n"
"} bind def\n"
"% draw bases\n"
"/drawbases {\n"
"  [] 0 setdash\n"
"  seqcolor\n"
"  0\n"
"  coor {\n"
"    aload pop moveto\n"
"    dup sequence exch 1 getinterval cshow\n"
"    1 add\n"
"  } forall\n"
"  pop\n"
"} bind def\n\n"
"/init {\n"
"  /Helvetica findfont fsize scalefont setfont\n"
"  1 setlinejoin\n"
"  1 setlinecap\n"
"  0.8 setlinewidth\n"
"  72 216 translate\n"
"  % find the coordinate range\n"
"  /xmax -1000 def /xmin 10000 def\n"
"  /ymax -1000 def /ymin 10000 def\n"
"  coor {\n"
"      aload pop\n"
"      dup ymin lt {dup /ymin exch def} if\n"
"      dup ymax gt {/ymax exch def} {pop} ifelse\n"
"      dup xmin lt {dup /xmin exch def} if\n"
"      dup xmax gt {/xmax exch def} {pop} ifelse\n"
"  } forall\n"
"  /size {xmax xmin sub ymax ymin sub max} bind def\n"
"  72 6 mul size div dup scale\n"
"  size xmin sub xmax sub 2 div size ymin sub ymax sub 2 div\n"
"  translate\n"
"} bind def\n"
"end\n";

static const char *RNAss_color =
"/range 0.8 def\n"
"/drawreliability {\n"
"  /Smax 1 def\n"
"  0\n"
"  coor {\n"
"    aload pop\n"
"    S 3 index get\n"
"    Smax div range mul\n"
"    invert {range exch sub} if\n"
"    1 1 sethsbcolor\n"
"    newpath\n"
"    fsize 2 div 0 360 arc\n"
"    fill\n"
"    1 add\n"
"  } forall\n"
"} bind def\n"
"/colorbar { %% xloc yloc colorbar -> []\n"
"  /STR 8 string def\n"
"  gsave\n"
"    xmin xmax add size sub 2 div\n"
"    ymin ymax add size sub 2 div translate\n"
"    size dup scale\n"
"    translate\n"
"    0.015 dup scale\n"
"    /tics 64 def\n"
"    gsave\n"
"      10 tics div 1 scale\n"
"      0 1 tics\n"
"      {\n"
"	  dup 0 moveto 0.5 add\n"
"	  tics div range mul\n"
"	  invert {range exch sub} if\n"
"	  1 1 sethsbcolor\n"
"	  1 0 rlineto 0 1 rlineto -1 0 rlineto closepath fill\n"
"      } for\n"
"    grestore\n"
"    0 setgray\n"
"    -0.1 1.01 moveto (0) gsave 0.1 dup scale show grestore\n"
"    10 1.01 moveto Smax STR cvs\n"
"    gsave 0.1 dup scale dup stringwidth pop -2 div 0 rmoveto show grestore\n"
"  grestore\n"
"} bind def\n";

// from Vienna RNA package 1.8.3
//int PS_rna_plot_a(char *string, char *structure, char *ssfile, char *pre, char *post)
int ps_plot(const char *string, const char *structure, const char* ssfile)
{
  float  xmin, xmax, ymin, ymax, size;
  int    i, length;
  float *X, *Y;
  FILE  *xyplot;
  short *pair_table;

  length = strlen(string);

  xyplot = fopen(ssfile, "w");
  if (xyplot == NULL) {
    fprintf(stderr, "can't open file %s - not doing xy_plot\n", ssfile);
    return 0;
  }

  pair_table = make_pair_table(structure);

  X = (float *) malloc((length+1)*sizeof(float));
  Y = (float *) malloc((length+1)*sizeof(float));
  i = naview_xy_coordinates(pair_table, X, Y);

  xmin = xmax = X[0];
  ymin = ymax = Y[0];
  for (i = 1; i < length; i++) {
     xmin = X[i] < xmin ? X[i] : xmin;
     xmax = X[i] > xmax ? X[i] : xmax;
     ymin = Y[i] < ymin ? Y[i] : ymin;
     ymax = Y[i] > ymax ? Y[i] : ymax;
  }
  size = MAX((xmax-xmin),(ymax-ymin));

  fprintf(xyplot,
	  "%%!PS-Adobe-3.0 EPSF-3.0\n"
	  "%%%%Creator: %s, ViennaRNA-%s\n"
	  "%%%%CreationDate: %s"
	  "%%%%Title: RNA Secondary Structure Plot\n"
	  "%%%%BoundingBox: 66 210 518 662\n"
	  "%%%%DocumentFonts: Helvetica\n"
	  "%%%%Pages: 1\n"
	  "%%%%EndComments\n\n"
	  "%%Options: %s\n", rcsid+5, VERSION, time_stamp(), "" /*option_string()*/);
  fprintf(xyplot, "%% to switch off outline pairs of sequence comment or\n"
	  "%% delete the appropriate line near the end of the file\n\n");
  fprintf(xyplot, "%s", RNAss_head);

#if 0
  if (pre || post) {
    fprintf(xyplot, "%s", anote_macros);
  }
#endif
  fprintf(xyplot, "%%%%EndProlog\n");

  fprintf(xyplot, "RNAplot begin\n"
	  "%% data start here\n");

#if 0
  /* cut_point */
  if (cut_point > 0 && cut_point <= strlen(string))
    fprintf(xyplot, "/cutpoint %d def\n", cut_point-1);
#endif
  /* sequence */
  fprintf(xyplot,"/sequence (\\\n");
  i=0;
  while (i<length) {
    fprintf(xyplot, "%.255s\\\n", string+i);  /* no lines longer than 255 */
    i+=255;
  }
  fprintf(xyplot,") def\n");
  /* coordinates */
  fprintf(xyplot, "/coor [\n");
  for (i = 0; i < length; i++)
    fprintf(xyplot, "[%3.3f %3.3f]\n", X[i], Y[i]);
  fprintf(xyplot, "] def\n");
  /* base pairs */
  fprintf(xyplot, "/pairs [\n");
  for (i = 1; i <= length; i++)
    if (pair_table[i]>i)
      fprintf(xyplot, "[%d %d]\n", i, pair_table[i]);
  fprintf(xyplot, "] def\n\n");

  fprintf(xyplot, "init\n\n");
  /* draw the data */
#if 0
  if (pre) {
    fprintf(xyplot, "%% Start Annotations\n");
    fprintf(xyplot, "%s\n", pre);
    fprintf(xyplot, "%% End Annotations\n");
  }
#endif
  fprintf(xyplot,
	  "%% switch off outline pairs or bases by removing these lines\n"
	  "drawoutline\n"
	  "drawpairs\n"
	  "drawbases\n");
#if 0
  if (post) {
    fprintf(xyplot, "%% Start Annotations\n");
    fprintf(xyplot, "%s\n", post);
    fprintf(xyplot, "%% End Annotations\n");
  }
#endif
  fprintf(xyplot, "%% show it\nshowpage\n");
  fprintf(xyplot, "end\n");
  fprintf(xyplot, "%%%%EOF\n");

  fclose(xyplot);

  free(pair_table);
  free(X); free(Y);
  return 1; /* success */
}

// from Vienna RNA package 1.8.3
//int PS_rna_plot_a(char *string, char *structure, char *ssfile, char *pre, char *post)
template < class BPTable >
int ps_color_plot(const char *string, const char *structure, const BPTable& bp, const char* ssfile)
{
  std::vector<typename BPTable::value_type> q(bp.calc_nonbp_prob());
  
  float  xmin, xmax, ymin, ymax, size;
  int    i, length;
  float *X, *Y;
  FILE  *xyplot;
  short *pair_table;

  length = strlen(string);

  xyplot = fopen(ssfile, "w");
  if (xyplot == NULL) {
    fprintf(stderr, "can't open file %s - not doing xy_plot\n", ssfile);
    return 0;
  }

  pair_table = make_pair_table(structure);

  X = (float *) malloc((length+1)*sizeof(float));
  Y = (float *) malloc((length+1)*sizeof(float));
  i = naview_xy_coordinates(pair_table, X, Y);

  xmin = xmax = X[0];
  ymin = ymax = Y[0];
  for (i = 1; i < length; i++) {
     xmin = X[i] < xmin ? X[i] : xmin;
     xmax = X[i] > xmax ? X[i] : xmax;
     ymin = Y[i] < ymin ? Y[i] : ymin;
     ymax = Y[i] > ymax ? Y[i] : ymax;
  }
  size = MAX((xmax-xmin),(ymax-ymin));

  fprintf(xyplot,
	  "%%!PS-Adobe-3.0 EPSF-3.0\n"
	  "%%%%Creator: %s, ViennaRNA-%s\n"
	  "%%%%CreationDate: %s"
	  "%%%%Title: RNA Secondary Structure Plot\n"
	  "%%%%BoundingBox: 66 210 518 662\n"
	  "%%%%DocumentFonts: Helvetica\n"
	  "%%%%Pages: 1\n"
	  "%%%%EndComments\n\n"
	  "%%Options: %s\n", rcsid+5, VERSION, time_stamp(), "" /*option_string()*/);
  fprintf(xyplot, "%% to switch off outline pairs of sequence comment or\n"
	  "%% delete the appropriate line near the end of the file\n\n");
  fprintf(xyplot, "%s", RNAss_head);

#if 0
  if (pre || post) {
    fprintf(xyplot, "%s", anote_macros);
  }
#endif
  fprintf(xyplot, "%%%%EndProlog\n");

  fprintf(xyplot, "RNAplot begin\n"
	  "%% data start here\n");

#if 0
  /* cut_point */
  if (cut_point > 0 && cut_point <= strlen(string))
    fprintf(xyplot, "/cutpoint %d def\n", cut_point-1);
#endif
  /* sequence */
  fprintf(xyplot,"/sequence (\\\n");
  i=0;
  while (i<length) {
    fprintf(xyplot, "%.255s\\\n", string+i);  /* no lines longer than 255 */
    i+=255;
  }
  fprintf(xyplot,") def\n");
  /* coordinates */
  fprintf(xyplot, "/coor [\n");
  for (i = 0; i < length; i++)
    fprintf(xyplot, "[%3.3f %3.3f]\n", X[i], Y[i]);
  fprintf(xyplot, "] def\n");
  /* base pairs */
  fprintf(xyplot, "/pairs [\n");
  for (i = 1; i <= length; i++)
    if (pair_table[i]>i)
      fprintf(xyplot, "[%d %d]\n", i, pair_table[i]);
  fprintf(xyplot, "] def\n\n");

  fprintf(xyplot, "init\n\n");
  /* draw the data */
#if 0
  if (pre) {
    fprintf(xyplot, "%% Start Annotations\n");
    fprintf(xyplot, "%s\n", pre);
    fprintf(xyplot, "%% End Annotations\n");
  }
#endif
  /* colored by base-pairing probabilities*/
  fprintf(xyplot, RNAss_color);
  fprintf(xyplot, "/S [\n");
  {
    std::vector<float> sp(length+1, -1.0);
    for (i = 1; i <= length; i++)
      if (pair_table[i]>i)
        sp[i] = sp[pair_table[i]] = bp(i-1, pair_table[i]-1);
    for (i = 1; i <= length; i++)
    {
      if (sp[i]<0.0) sp[i] = q[i-1];
      fprintf(xyplot, "  %7.5f\n", sp[i]);
    }
  }
  fprintf(xyplot, "] def\n\n");
  fprintf(xyplot, "/invert true def\n");
  fprintf(xyplot, "drawreliability\n");
  fprintf(xyplot, "0.0 0.0 colorbar\n");
  
  fprintf(xyplot,
	  "%% switch off outline pairs or bases by removing these lines\n"
	  "drawoutline\n"
	  "drawpairs\n"
	  "drawbases\n");
#if 0
  if (post) {
    fprintf(xyplot, "%% Start Annotations\n");
    fprintf(xyplot, "%s\n", post);
    fprintf(xyplot, "%% End Annotations\n");
  }
#endif
  fprintf(xyplot, "%% show it\nshowpage\n");
  fprintf(xyplot, "end\n");
  fprintf(xyplot, "%%%%EOF\n");

  fclose(xyplot);

  free(pair_table);
  free(X); free(Y);
  return 1; /* success */
}

// from Vienna RNA package 1.8.3
static short *make_pair_table(const char *structure)
{
    /* returns array representation of structure.
       table[i] is 0 if unpaired or j if (i.j) pair.  */
   short i,j,hx;
   short length;
   short *stack;
   short *table;

   length = (short) strlen(structure);
   stack = (short *) malloc(sizeof(short)*(length+1));
   table = (short *) malloc(sizeof(short)*(length+2));
   table[0] = length;

   for (hx=0, i=1; i<=length; i++) {
      switch (structure[i-1]) {
       case '(':
	 stack[hx++]=i;
	 break;
       case ')':
	 j = stack[--hx];
	 if (hx<0) {
	    fprintf(stderr, "%s\n", structure);
	    nrerror("unbalanced brackets in make_pair_table");
	 }
	 table[i]=j;
	 table[j]=i;
	 break;
       default:   /* unpaired base, usually '.' */
	 table[i]= 0;
	 break;
      }
   }
   if (hx!=0) {
      fprintf(stderr, "%s\n", structure);
      nrerror("unbalanced brackets in make_pair_table");
   }
   free(stack);
   return(table);
}

static char *time_stamp(void)
{
  time_t  cal_time;

  cal_time = time(NULL);
  return ( ctime(&cal_time) );
}

static void nrerror(const char message[])       /* output message upon error */
{
  fprintf(stderr, "\n%s\n", message);
  exit(EXIT_FAILURE);
}

#include "bp.h"
template
int ps_color_plot(const char *string, const char *structure, const BPTableTmpl<float>& bp, const char* ssfile);
