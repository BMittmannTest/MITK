/*===================================================================

The Medical Imaging Interaction Toolkit (MITK)

Copyright (c) German Cancer Research Center,
Division of Medical and Biological Informatics.
All rights reserved.

This software is distributed WITHOUT ANY WARRANTY; without
even the implied warranty of MERCHANTABILITY or FITNESS FOR
A PARTICULAR PURPOSE.

See LICENSE.txt or http://www.mitk.org for details.

===================================================================*/
// adapted from Nicolas Toussaint and Pierre Fillard vtkINRIA3D: A VTK Extension for
// Spatiotemporal Data Synchronization, Visualization and Management. (CECILL-2.0)
#ifndef _lut_BlackBody_h_
#define _lut_BlackBody_h_


static const int BlackBody[256*3]={
0,
3,
6,
9,
12,
15,
18,
21,
24,
27,
30,
33,
36,
39,
42,
45,
48,
51,
54,
57,
60,
63,
66,
69,
72,
75,
78,
81,
84,
87,
90,
93,
96,
99,
102,
105,
108,
111,
114,
117,
120,
123,
126,
129,
132,
135,
138,
141,
144,
147,
150,
153,
156,
159,
162,
165,
168,
171,
174,
177,
180,
183,
186,
189,
192,
195,
198,
201,
204,
207,
210,
213,
216,
219,
222,
225,
228,
231,
234,
237,
240,
243,
246,
249,
252,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
3,
6,
9,
12,
15,
18,
21,
24,
27,
30,
33,
36,
39,
42,
45,
48,
51,
54,
57,
60,
63,
66,
69,
72,
75,
78,
81,
84,
87,
90,
93,
96,
99,
102,
105,
108,
111,
114,
117,
120,
123,
126,
129,
132,
135,
138,
141,
144,
147,
150,
153,
156,
159,
162,
165,
168,
171,
174,
177,
180,
183,
186,
189,
192,
195,
198,
201,
204,
207,
210,
213,
216,
219,
222,
225,
228,
231,
234,
237,
240,
243,
246,
249,
252,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
255,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
0,
3,
6,
9,
12,
15,
18,
21,
24,
27,
30,
33,
36,
39,
42,
45,
48,
51,
54,
57,
60,
63,
66,
69,
72,
75,
78,
81,
84,
87,
90,
93,
96,
99,
102,
105,
108,
111,
114,
117,
120,
123,
126,
129,
132,
135,
138,
141,
144,
147,
150,
153,
156,
159,
162,
165,
168,
171,
174,
177,
180,
183,
186,
189,
192,
195,
198,
201,
204,
207,
210,
213,
216,
219,
222,
225,
228,
231,
234,
237,
240,
243,
246,
249,
252,
255
};

#endif
