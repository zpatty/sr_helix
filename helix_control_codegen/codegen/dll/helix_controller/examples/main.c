/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 * File: main.c
 *
 * MATLAB Coder version            : 23.2
 * C/C++ source code generated on  : 28-Oct-2024 17:40:15
 */

/*************************************************************************/
/* This automatically generated example C main file shows how to call    */
/* entry-point functions that MATLAB Coder generated. You must customize */
/* this file for your application. Do not modify this file directly.     */
/* Instead, make a copy of this file, modify it, and integrate it into   */
/* your development environment.                                         */
/*                                                                       */
/* This file initializes entry-point function arguments to a default     */
/* size and value before calling the entry-point functions. It does      */
/* not store or use any values returned from the entry-point functions.  */
/* If necessary, it does pre-allocate memory for returned values.        */
/* You can use this file as a starting point for a main function that    */
/* you can deploy in your application.                                   */
/*                                                                       */
/* After you copy the file, and before you deploy it, you must make the  */
/* following changes:                                                    */
/* * For variable-size function arguments, change the example sizes to   */
/* the sizes that your application requires.                             */
/* * Change the example values of function arguments to the values that  */
/* your application requires.                                            */
/* * If the entry-point functions return values, store these values or   */
/* otherwise use them as required by your application.                   */
/*                                                                       */
/*************************************************************************/

/* Include Files */
#include "main.h"
#include "helix_controller.h"
#include "helix_controller_emxAPI.h"
#include "helix_controller_terminate.h"
#include "helix_controller_types.h"
#include "rt_nonfinite.h"

/* Function Declarations */
static void argInit_10x10_real_T(double result[100]);

static void argInit_10x1_real_T(double result[10]);

static void argInit_3x1_real_T(double result[3]);

static double argInit_real_T(void);

/* Function Definitions */
/*
 * Arguments    : double result[100]
 * Return Type  : void
 */
static void argInit_10x10_real_T(double result[100])
{
  int i;
  /* Loop over the array to initialize each element. */
  for (i = 0; i < 100; i++) {
    /* Set the value of the array element.
Change this value to the value that the application requires. */
    result[i] = argInit_real_T();
  }
}

/*
 * Arguments    : double result[10]
 * Return Type  : void
 */
static void argInit_10x1_real_T(double result[10])
{
  int idx0;
  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 10; idx0++) {
    /* Set the value of the array element.
Change this value to the value that the application requires. */
    result[idx0] = argInit_real_T();
  }
}

/*
 * Arguments    : double result[3]
 * Return Type  : void
 */
static void argInit_3x1_real_T(double result[3])
{
  int idx0;
  /* Loop over the array to initialize each element. */
  for (idx0 = 0; idx0 < 3; idx0++) {
    /* Set the value of the array element.
Change this value to the value that the application requires. */
    result[idx0] = argInit_real_T();
  }
}

/*
 * Arguments    : void
 * Return Type  : double
 */
static double argInit_real_T(void)
{
  return 0.0;
}

/*
 * Arguments    : int argc
 *                char **argv
 * Return Type  : int
 */
int main(int argc, char **argv)
{
  (void)argc;
  (void)argv;
  /* The initialize function is being called automatically from your entry-point
   * function. So, a call to initialize is not included here. */
  /* Invoke the entry-point functions.
You can call entry-point functions multiple times. */
  main_helix_controller();
  /* Terminate the application.
You do not need to do this more than one time. */
  helix_controller_terminate();
  return 0;
}

/*
 * Arguments    : void
 * Return Type  : void
 */
void main_helix_controller(void)
{
  emxArray_real_T *A;
  emxArray_real_T *C;
  emxArray_real_T *M;
  emxArray_real_T *cq;
  emxArray_real_T *tau;
  double Kp_tmp[100];
  double x[36];
  double q_tmp[10];
  double xd_tmp[3];
  double N_tmp;
  double tau_r;
  /* Initialize function 'helix_controller' input arguments. */
  /* Initialize function input argument 'q'. */
  argInit_10x1_real_T(q_tmp);
  /* Initialize function input argument 'dq'. */
  /* Initialize function input argument 'qd'. */
  /* Initialize function input argument 'dqd'. */
  /* Initialize function input argument 'ddqd'. */
  N_tmp = argInit_real_T();
  /* Initialize function input argument 'Kp'. */
  argInit_10x10_real_T(Kp_tmp);
  /* Initialize function input argument 'KD'. */
  /* Initialize function input argument 'xd'. */
  argInit_3x1_real_T(xd_tmp);
  /* Initialize function input argument 'dxd'. */
  /* Initialize function input argument 'dxr'. */
  /* Call the entry-point 'helix_controller'. */
  emxInitArray_real_T(&tau, 1);
  emxInitArray_real_T(&M, 2);
  emxInitArray_real_T(&C, 1);
  emxInitArray_real_T(&A, 2);
  emxInitArray_real_T(&cq, 1);
  helix_controller(q_tmp, q_tmp, q_tmp, q_tmp, q_tmp, N_tmp, N_tmp, N_tmp,
                   N_tmp, N_tmp, N_tmp, N_tmp, N_tmp, N_tmp, N_tmp, Kp_tmp,
                   Kp_tmp, N_tmp, N_tmp, xd_tmp, xd_tmp, xd_tmp, N_tmp, N_tmp,
                   tau, &tau_r, x, M, C, A, cq);
  emxDestroyArray_real_T(tau);
  emxDestroyArray_real_T(M);
  emxDestroyArray_real_T(C);
  emxDestroyArray_real_T(A);
  emxDestroyArray_real_T(cq);
}

/*
 * File trailer for main.c
 *
 * [EOF]
 */
