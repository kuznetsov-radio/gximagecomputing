#pragma once

#define D2(s1, i1, i2)                 ((i1)+(i2)*(s1))
#define D3(s1, s2, i1, i2, i3)         ((i1)+((i2)+(i3)*(s2))*(s1))
#define D4(s1, s2, s3, i1, i2, i3, i4) ((i1)+((i2)+((i3)+(i4)*(s3))*(s2))*(s1))