/*
 * lamb.h — Gooding Lambert solver
 * Original by R.H. Gooding, adapted by N. Strange
 */

#ifndef LAMB_H
#define LAMB_H

int lambert(double gm, double r1[3], double r2[3], int nrev, double dt,
            double v1[3], double v2[3]);

#endif
