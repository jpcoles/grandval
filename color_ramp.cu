#include <math.h>
#include "color_ramp.h"

void color_ramp_grey(float *r, float *g, float *b)
{
    const float r0 = *r;
    const float g0 = *g;
    const float b0 = *b;

    float c = sqrtf(r0*r0 + g0*g0 + b0*b0);
    if (c > 1.0F) c = 1.0F;

    *r = *g = *b = c;
}

void color_ramp_hot2cold(float *r, float *g, float *b)
{
    const float r0 = *r;
    const float g0 = *g;
    const float b0 = *b;

    float cr, cg, cb;
    cr = cg = cb = 1.0F;

    float v = sqrtf(r0*r0 + g0*g0 + b0*b0);

    if (v < 0.0F)
      v = 0;
    if (v > 1.0F)
      v = 1.0F;

    if (v < 0.25) {
      cr = 0;
      cg = 4 * v;
    } else if (v < (0.5)) {
      cr = 0;
      cb = 1 + 4 * (0.25 - v);
    } else if (v < (0.75)) {
      cr = 4 * (v - 0.5);
      cb = 0;
    } else {
      cg = 1 + 4 * (0.75 - v);
      cb = 0;
    }

    *r = cr;
    *g = cg;
    *b = cb;

}

void color_ramp_astro(float *r, float *g, float *b)
{
    const float r0 = *r;
    const float g0 = *g;
    const float b0 = *b;
    float v = sqrtf(r0*r0 + g0*g0 + b0*b0);

    static int ramp[256][3] =  {
  {0,0,255}, {0,0,254}, {1,0,253}, {2,0,252}, {3,0,250}, {4,0,249}, {5,0,248}, {6,0,247}, {7,0,247},
  {8,0,246}, {9,0,245}, {10,0,244}, {11,0,242}, {12,0,241}, {13,0,240}, {14,0,239}, {15,0,239}, {16,0,238},
 {17,0,237}, {18,0,236}, {19,0,234}, {20,0,233}, {21,0,232}, {22,0,231}, {23,0,231}, {24,0,230}, {25,0,229},
 {26,0,228}, {27,0,226}, {28,0,225}, {29,0,224}, {30,0,223}, {31,0,223}, {33,0,222}, {33,0,221}, {35,0,220},
 {35,0,218}, {37,0,217}, {37,0,216}, {39,0,215}, {39,0,215}, {41,0,214}, {41,0,213}, {43,0,212}, {43,0,210},
 {45,0,209}, {45,0,208}, {47,0,207}, {47,0,207}, {49,0,206}, {49,0,205}, {51,0,204}, {51,0,202}, {53,0,201},
 {53,0,200}, {55,0,199}, {55,0,199}, {57,0,198}, {57,0,197}, {59,0,196}, {59,0,194}, {61,0,193}, {61,0,192},
 {63,0,191}, {63,0,191}, {64,0,190}, {66,0,189}, {67,0,188}, {67,0,186}, {68,0,185}, {70,0,184}, {71,0,183},
 {71,0,183}, {72,0,182}, {74,0,181}, {75,0,180}, {75,0,178}, {76,0,177}, {78,0,176}, {79,0,175}, {79,0,175},
 {80,0,174}, {82,0,173}, {83,0,172}, {83,0,170}, {84,0,169}, {86,0,168}, {87,0,167}, {87,0,167}, {88,0,166},
 {90,0,165}, {91,0,164}, {91,0,162}, {92,0,161}, {94,0,160}, {95,0,159}, {95,0,159}, {96,0,158}, {98,0,157},
 {99,0,156}, {99,0,154}, {100,0,153}, {102,0,152}, {103,0,151}, {103,0,151}, {104,0,150}, {106,0,149}, {107,0,148},
{107,0,146}, {108,0,145}, {110,0,144}, {111,0,143}, {111,0,143}, {112,0,142}, {114,0,141}, {115,0,140}, {115,0,138},
{116,0,137}, {118,0,136}, {119,0,135}, {119,0,135}, {120,0,134}, {122,0,133}, {123,0,132}, {123,0,130}, {124,0,129},
{126,0,128}, {127,0,127}, {127,0,127}, {128,0,126}, {129,0,124}, {130,0,123}, {132,0,123}, {133,0,122}, {134,0,120},
{135,0,119}, {135,0,119}, {136,0,118}, {137,0,116}, {138,0,115}, {140,0,115}, {141,0,114}, {142,0,112}, {143,0,111},
{143,0,111}, {144,0,110}, {145,0,108}, {146,0,107}, {148,0,107}, {149,0,106}, {150,0,104}, {151,0,103}, {151,0,103},
{152,0,102}, {153,0,100}, {154,0, 99}, {156,0, 99}, {157,0, 98}, {158,0, 96}, {159,0, 95}, {159,0, 95}, {160,0, 94},
{161,0, 92}, {162,0, 91}, {164,0, 91}, {165,0, 90}, {166,0, 88}, {167,0, 87}, {167,0, 87}, {168,0, 86}, {169,0, 84},
{170,0, 83}, {172,0, 83}, {173,0, 82}, {174,0, 80}, {175,0, 79}, {175,0, 79}, {176,0, 78}, {177,0, 76}, {178,0, 75},
{180,0, 75}, {181,0, 74}, {182,0, 72}, {183,0, 71}, {183,0, 71}, {184,0, 70}, {185,0, 68}, {186,0, 67}, {188,0, 67},
{189,0, 66}, {190,0, 64}, {191,0, 63}, {191,0, 63}, {192,0, 61}, {193,0, 61}, {194,0, 59}, {196,0, 59}, {197,0, 57},
{198,0, 57}, {199,0, 55}, {199,0, 55}, {200,0, 53}, {201,0, 53}, {202,0, 51}, {204,0, 51}, {205,0, 49}, {206,0, 49},
{207,0, 47}, {207,0, 47}, {208,0, 45}, {209,0, 45}, {210,0, 43}, {212,0, 43}, {213,0, 41}, {214,0, 41}, {215,0, 39},
{215,0, 39}, {216,0, 37}, {217,0, 37}, {218,0, 35}, {220,0, 35}, {221,0, 33}, {222,0, 33}, {223,0, 31}, {223,0, 30},
{224,0, 29}, {225,0, 28}, {226,0, 27}, {228,0, 26}, {229,0, 25}, {230,0, 24}, {231,0, 23}, {231,0, 22}, {232,0, 21},
{233,0, 20}, {234,0, 19}, {236,0, 18}, {237,0, 17}, {238,0, 16}, {239,0, 15}, {239,0, 14}, {240,0, 13}, {241,0, 12},
{242,0, 11}, {244,0, 10}, {245,0,  9}, {246,0,  8}, {247,0,  7}, {247,0,  6}, {248,0,  5}, {249,0,  4}, {250,0,  3},
{252,0,  2}, {253,0,  1}, {254,0,  0}, {255,0,  0}};

    int i = v * 256;
    if (i > 255) i = 255;
    *r = ramp[i][0] / 255.0F;
    *g = ramp[i][1] / 255.0F;
    *b = ramp[i][2] / 255.0F;
}

void color_ramp_wrbb(float *r, float *g, float *b)
{
    static int first_time = 1;
    static double ramp[256][3];
    double slope, offset;
    int i,j;

    if (first_time)
    {
        first_time = 0;
        slope = 255./20.;
        for(i = 0 ;i < 21 ;i++){
            ramp[i][0] = 0 ;
            ramp[i][1] = 0 ;
            ramp[i][2] = (int)(slope * (double)(i) + .5) ;
        }
        slope = 191./20. ;
        offset = 64. - slope * 21. ;
        for(i = 21 ;i < 42 ;i++){
            ramp[i][0] = (int)(slope * (double)(i) + offset + .5) ;
            ramp[i][1] = 0 ;
            ramp[i][2] = 255 ;
        }
        slope = -205./21. ;
        offset = 255. - slope * 41. ;
        for(i = 42 ;i < 63 ;i++){
            ramp[i][0] = 255 ;
            ramp[i][1] = 0 ;
            ramp[i][2] = (int)(slope * (double)(i) + offset + .5) ;
        }
        slope = 205./40. ;
        offset = 50. - slope * 63. ;
        for(i = 63 ;i < 104 ;i++){
            ramp[i][0] = 255 ;
            ramp[i][1] = (int)(slope * (double)(i) + offset + .5) ;
            ramp[i][2] = 0 ;
        }
        slope = 255./21. ;
        offset = -slope * 103. ;
        for(i = 104 ;i < 125 ;i++){
            ramp[i][0] = 255 ;
            ramp[i][1] = 255 ;
            ramp[i][2] = (int)(slope * (double)(i) + offset + .5) ;
        }
        for (i = 125, j = 0 ; j < 125 ;i++, j++) {
            ramp[i][0] = ramp[j][0];
            ramp[i][1] = ramp[j][1];
            ramp[i][2] = ramp[j][2];
        }
    }

    const float r0 = *r;
    const float g0 = *g;
    const float b0 = *b;
    float v = sqrtf(r0*r0 + g0*g0 + b0*b0);
    i = v * 256;
    if (i > 255) i = 255;
    *r = ramp[i][0] / 255.0F;
    *g = ramp[i][1] / 255.0F;
    *b = ramp[i][2] / 255.0F;
}
