/* ************************************************************************
 * Copyright 2013 Advanced Micro Devices, Inc.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 * http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 * ************************************************************************/


#include <sys/types.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

/* Include CLBLAS header. It automatically includes needed OpenCL header,
 * so we can drop out explicit inclusion of cl.h header.
 */
#include <clBLAS.h>

//#include "dspBLAS_stub.h"

/* This example uses predefined matrices and their characteristics for
 * simplicity purpose.
 */

#define M  4
#define N  3
#define K  5

static const clblasOrder order = clblasRowMajor;

static const cl_float alpha = 1;

static const clblasTranspose transA = clblasNoTrans;
static const cl_float A[M*K] = {
    11, 12, 13, 14, 15,
    21, 22, 23, 24, 25,
    31, 32, 33, 34, 35,
    41, 42, 43, 44, 45,
};
static const size_t lda = K;        /* i.e. lda = K */

static const clblasTranspose transB = clblasNoTrans;
static const cl_float B[K*N] = {
    11, 12, 13,
    21, 22, 23,
    31, 32, 33,
    41, 42, 43,
    51, 52, 53,
};
static const size_t ldb = N;        /* i.e. ldb = N */

static const cl_float beta = 1;

static cl_float C[M*N] = {
    11, 12, 13,
    21, 22, 23,
    31, 32, 33,
    41, 42, 43,
};
static const size_t ldc = N;        /* i.e. ldc = N */

static cl_float result[M*N];

static const size_t off  = 0;
static const size_t offA = 0;//K + 1;   /* K + off */
static const size_t offB = 0;//N + 1;   /* N + off */
static const size_t offC = 0;//N + 1;   /* N + off */

static void
printResult(const char* str)
{
    size_t i, j, nrows;

    printf("%s:\n", str);

int ldc_ = ldc;
//int ldc_ = K;
    nrows = (sizeof(result) / sizeof(cl_float)) / ldc_;
    for (i = 0; i < nrows; i++) {
        for (j = 0; j < ldc_; j++) {
            printf("%d ", (int)result[i * ldc_ + j]);
        }
        printf("\n");
    }
}

int m, k, n;
#if 0
static void
printMatrixf32(const char *msg, float *m, size_t len, size_t ld)
{
    size_t i, j, nrows;
    nrows = len / ld;
printf("nrows=%u len=%u ld=%u\n", nrows, len, ld);
    printf("<< %s >>\n", msg);
    for (i = 0; i < nrows; i++) {
        for (j = 0; j < ld; j++) {
            printf("%d ", (int)m[i * ld + j]);
        }
        printf("\n");
    }
}

int do_dsp_gemm(int m, int n, int k)
{
	int i;
	int loop = 20;
	int ret = 0;
	float *a = NULL, *b = NULL, *c = NULL;
	float alpha = 1, beta = 1;

	hexagon_init();

	a = (float *)dspmem_alloc(NULL, 0, m * k * sizeof(float), 16);
	b = (float *)dspmem_alloc(NULL, 0, n * k * sizeof(float), 16);
	c = (float *)dspmem_alloc(NULL, 0, m * n * sizeof(float), 16);
	if (!a || !b || !c) {
		ret = -1;
		goto exit;
	}
	for (i = 0; i < m * k; ++i) {
		a[i] = i % 100;
	}
	for (i = 0; i < n * k; ++i) {
		b[i] = i % 50;
	}
	for (i = 0; i < m * n; ++i) {
		c[i] = i % 20;
	}
	for (i = 0; i < loop; ++i) {
		dspBLAS_SGEMMf32(0, 0, 0,
				 m, n, k,
				 alpha,
				 a, m * k, 0, k,
				 b, n * k, 0, n,
				 beta,
				 c, m * k,
				 c, m * k, 0, n);
	}
exit:
	if (a) {
		dspmem_free(a);
	}
	if (b) {
		dspmem_free(b);
	}
	if (c) {
		dspmem_free(c);
	}
	hexagon_deinit();
	return ret;
}

#define LEN(_m_)  (sizeof(_m_) / sizeof(_m_[0]))
void
call_dspBLAS(void)
{
	float *bufA, *bufB, *bufC;
	float *bufA_;
	int ret;
	int pivot = 2;


do_dsp_gemm(m, n, k);
return;
	//int m = M;
	//int k = K;
	//int n = N;
	//int m = 50;
	//int k = 500;
	//int n = 64;

	hexagon_init();
	dspBLAS_SanityTest(&pivot);
	printf("pivot = %d\n", pivot);
#if 1
printf("%s:%d\n", __FUNCTION__, __LINE__);
	bufA = (float *)dspmem_alloc(NULL, 0, sizeof(float) * m * k, 16);
printf("%s:%d\n", __FUNCTION__, __LINE__);
	bufA_ = (float *)dspmem_alloc(NULL, 0, sizeof(float) * m * k, 16);
	bufB = (float *)dspmem_alloc(NULL, 0, sizeof(float) * k * n, 16);
	bufC = (float *)dspmem_alloc(NULL, 0, sizeof(float) * m * n * k, 16);
	//bufA = (float *)dspmem_alloc(NULL, 0, sizeof(A), 16);
	//bufA_ = (float *)dspmem_alloc(NULL, 0, sizeof(A), 16);
	//bufB = (float *)dspmem_alloc(NULL, 0, sizeof(B), 16);
	//bufC = (float *)dspmem_alloc(NULL, 0, sizeof(C), 16);
	if (!bufA || !bufB || !bufC) {
		goto error;
	}
	printf("bufA = %p, bufA_ = %p\n", bufA, bufA_);
	printf("bufB = %p, bufC = %p\n", bufB, bufC);

	// 128-bit aligned addresses
	assert(!((uint32_t)bufA & 0x0f));
	assert(!((uint32_t)bufB & 0x0f));
	assert(!((uint32_t)bufC & 0x0f));
	assert(!((uint32_t)bufA_ & 0x0f));

	memcpy(bufA, A, sizeof(A));
	memcpy(bufB, B, sizeof(B));
	memcpy(bufC, C, sizeof(C));
#else
	bufA = (float *)A;
	bufB = (float *)B;
	bufC = C;
#endif
	memset(result, 0, sizeof(C));
	printMatrixf32("bufA", bufA, LEN(A), lda);
	printMatrixf32("bufB", bufB, LEN(B), ldb);
	printMatrixf32("bufC", bufC, LEN(C), ldc);
#if 1
	int ldc_ = n;
	int m_ = m, n_ = n, k_ = k;


	#define min(a, b) (a > b ? (b) : (a))
	printf("=== %d %d %d\n", m, n, k);
	//printf("Memory  io: %d\n", (m * k) * 4 * 2);
	//if (transA == clblasTrans) {
	//	int ret = DSP_VERIFY(dspBLAS_Transposef32(bufA, m * k, k, m, bufA_, m * k));
	//	if (ret < 0) {
	//		goto error;
	//	}
	//	int tmp = m;
	//	m = k;
	//	k = min(k, tmp);
	//}
	//if (transB == clblasTrans) {
	//	int tmp = n;
	//	n = K;
	//	k = min(k, tmp);

	//}
	//printMatrixf32("bufA_", bufA_, LEN(A), m);
	//printf("Memory  in: %6d\n", (m * k + k * n) * 4);
	//printf("Memory out: %6d\n", 4 * (m * n));
	//printf("Memory tot: %6d\n", 4 * ((m * n) + m * k + k * n));
	//printf("after trans ===> %d %d %d\n", m, n, k);

	// row major
	int ta, tb;
	ta = transA == clblasTrans ? (1) : (0);
	tb = transB == clblasTrans ? (1) : (0),
	printf("transA = %d, transB = %d\n", ta, tb);
	for (int i = 0; i < 1; ++i) {
		int ret = 0;
		printf("iter=%d\n", i);
		//DSP_VERIFY(dspBLAS_SGEMMf32(order, ta, tb,
		//	         m, n, k,
		//		 alpha,
		//		 bufA, LEN(A), offA, lda, bufB, LEN(B), offB, ldb,
		//		 beta,
		//		 bufC, LEN(C), offC, ldc));
#if 1
		DSP_VERIFY(dspBLAS_SGEMMf32(0, ta, tb,
				 m, n, k,
				 alpha,
				 bufA, m * k, 0, k,
				 bufB, k * n, 0, n,
				 beta,
				 bufC, m * n,
				 bufC, m * n, 0, n), &ret);
		DSP_VERIFY(dspBLAS_repeat(), &ret);
#else
		int am = m, ak = k, as = k;
		int bn = n, bk = k, bs = n;
		int cm = m, cn = n;
		if (alpha != 0) {
			if (transA == clblasTrans) {
				DSP_VERIFY(dspBLAS_Transposef32(
					bufA, m * k, k, m, 0, bufA, m * k, 0), 0);
				printMatrixf32("T(bufA)", bufA, LEN(A), M);
				am = k;
				ak = min(m, k);
				as = m;
				cm = k;
			}
			if (transB == clblasTrans) {
				DSP_VERIFY(dspBLAS_Transposef32(
					bufB, k * n, n, k, 0, bufB, k * n, 0), 0);
				bn = k;
				bk = n;
				ak = min(bk, ak);
				printMatrixf32("T(bufB)", bufB, LEN(B), K);
				cn = k;
			}
			DSP_VERIFY(dspBLAS_MatrixMultiplyf32(
					 bufA, m * k, ak, am, as * sizeof(float),
					 bufB, k * n, bn, bs * sizeof(float),
					 bufC, cm * cn, 0), &ret);
			if (alpha != 1) {
				DSP_VERIFY(dspBLAS_MultiplyScalarf32(
					 bufC, cm * cn, cn, cm, 0, alpha, bufC, cm * cn, 0), 0);
			}
		}
		if (beta != 0) {
			float *bufC_ = dspmem_alloc(0, 0, n * m * sizeof(float), 16);
			memcpy(bufC_, C, sizeof(C));
			DSP_VERIFY(dspBLAS_MultiplyScalarf32(
				 bufC_, m * n, n, m, 0, beta, bufC_, m * n, 0), 0);
			printMatrixf32("bufC_", bufC_, LEN(C), n);
			DSP_VERIFY(dspBLAS_Addf32(
				 bufC, m * n, n, m, 0, bufC_, m * n, 0, bufC, m * n, 0), 0);
			dspmem_free(bufC_);
		}
#endif
		if (ret < 0) {
			goto error;
		}
	 }
	if (transA == clblasTrans) {
		int tmp = m;
		m = k;
		k = min(k, tmp);
	}
	if (transB == clblasTrans) {
		int tmp = n;
		n = K;
		k = min(k, tmp);
		ldc_ = K;
	}
#endif
	printMatrixf32("bufC", bufC, ldc_ * m, ldc_);
	pivot = 2;
	dspBLAS_SanityTest(&pivot);
	printf("pivot = %d\n", pivot);
error:
	dspmem_free((void *)bufA);
	dspmem_free((void *)bufA_);
	dspmem_free((void *)bufB);
	dspmem_free((void *)bufC);
	hexagon_deinit();
}
#endif

int
main(int argc, char *argv[])
{
    cl_int err;
    cl_platform_id platform = 0;
    cl_device_id device = 0;
    cl_context_properties props[3] = { CL_CONTEXT_PLATFORM, 0, 0 };
    cl_context ctx = 0;
    cl_command_queue queue = 0;
    cl_mem bufA, bufB, bufC;
    cl_event event = NULL;
    int ret = 0;

    /* Setup OpenCL environment. */
    err = clGetPlatformIDs(1, &platform, NULL);
    if (err != CL_SUCCESS) {
        printf( "clGetPlatformIDs() failed with %d\n", err );
        return 1;
    }

    err = clGetDeviceIDs(platform, CL_DEVICE_TYPE_GPU, 1, &device, NULL);
    if (err != CL_SUCCESS) {
        printf( "clGetDeviceIDs() failed with %d\n", err );
        return 1;
    }

    props[1] = (cl_context_properties)platform;
    ctx = clCreateContext(props, 1, &device, NULL, NULL, &err);
    if (err != CL_SUCCESS) {
        printf( "clCreateContext() failed with %d\n", err );
        return 1;
    }

    queue = clCreateCommandQueue(ctx, device, 0, &err);
    if (err != CL_SUCCESS) {
        printf( "clCreateCommandQueue() failed with %d\n", err );
        clReleaseContext(ctx);
        return 1;
    }

    /* Setup clblas. */
    err = clblasSetup();
    if (err != CL_SUCCESS) {
        printf("clblasSetup() failed with %d\n", err);
        clReleaseCommandQueue(queue);
        clReleaseContext(ctx);
        return 1;
    }

    /* Prepare OpenCL memory objects and place matrices inside them. */
    bufA = clCreateBuffer(ctx, CL_MEM_READ_ONLY, M * K * sizeof(*A),
                          NULL, &err);
    bufB = clCreateBuffer(ctx, CL_MEM_READ_ONLY, K * N * sizeof(*B),
                          NULL, &err);
    bufC = clCreateBuffer(ctx, CL_MEM_READ_WRITE, M * N * sizeof(*C),
                          NULL, &err);

    err = clEnqueueWriteBuffer(queue, bufA, CL_TRUE, 0,
        M * K * sizeof(*A), A, 0, NULL, NULL);
    err = clEnqueueWriteBuffer(queue, bufB, CL_TRUE, 0,
        K * N * sizeof(*B), B, 0, NULL, NULL);
    err = clEnqueueWriteBuffer(queue, bufC, CL_TRUE, 0,
        M * N * sizeof(*C), C, 0, NULL, NULL);

    /* Call clblas extended function. Perform gemm for the lower right sub-matrices */
    err = clblasSgemm(order, transA, transB, M - off, N - off, K - off,
                         alpha, bufA, offA, lda,
                         bufB, offB, ldb, beta,
                         bufC, offC, ldc,
                         1, &queue, 0, NULL, &event);
    if (err != CL_SUCCESS) {
        printf("clblasSgemmEx() failed with %d\n", err);
        ret = 1;
    }
    else {
        /* Wait for calculations to be finished. */
        err = clWaitForEvents(1, &event);

        /* Fetch results of calculations from GPU memory. */
        err = clEnqueueReadBuffer(queue, bufC, CL_TRUE, 0,
                                  M * N * sizeof(*result),
                                  result, 0, NULL, NULL);

        /* At this point you will get the result of SGEMM placed in 'result' array. */
        puts("");
        printResult("clblasSgemmEx result");
    }

    /* Release OpenCL memory objects. */
    clReleaseMemObject(bufC);
    clReleaseMemObject(bufB);
    clReleaseMemObject(bufA);

    /* Finalize work with clblas. */
    clblasTeardown();

    /* Release OpenCL working objects. */
    clReleaseCommandQueue(queue);
    clReleaseContext(ctx);

    printf("=== %s %s %s\n", argv[1], argv[2], argv[3]);
    m = strtoul(argv[1], NULL, 10);
    n = strtoul(argv[2], NULL, 10);
    k = strtoul(argv[3], NULL, 10);
    printf("=== %d %d %d\n", m, n, k);

    // dsp
    //call_dspBLAS();


    return ret;
}
