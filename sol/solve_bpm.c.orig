#include "obpm.h"
#include "complex.h"
#include "obpm_prototype.h"
#include "ev.h"

#ifdef __NVCC__
  #include <thrust/complex.h>
  typedef thrust::complex<float> floatcomplex;
  #define I thrust::complex<float>{0,1}
  #define CEXPF(x) (thrust::exp(x))
  #define CREALF(x) (x.real())
  #define CIMAGF(x) (x.imag())
  #define MAX(x,y) (max(x,y))
  #define MIN(x,y) (min(x,y))
  #define FLOORF(x) (floor(x))
  #include <nvml.h>
  #define TILE_DIM 32
#else
  #ifdef __GNUC__ // This is defined for GCC and CLANG but not for Microsoft Visual C++ compiler
    //#define MAX(a,b) ({__typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a > _b? _a: _b;})
    //#define MIN(a,b) ({__typeof__ (a) _a = (a); __typeof__ (b) _b = (b); _a > _b? _b: _a;})
    // Prevent conflicting definition of d_complex_t
    #include <stdbool.h>
    #include <complex.h>

    //typedef float complex floatcomplex;
    typedef float _Complex floatcomplex;
    //#define I (1.0f * _Complex_I)
    // 手動で虚数単位を定義
    #define I ((float _Complex) 0.0f + 1.0f * (float _Complex) 1.0f)
    #define CEXPF(x) (cexpf(x))
    #define CREALF(x) (crealf(x))
    #define CIMAGF(x) (cimagf(x))
    //#define I std::complex<float>{0,1}
    //#define CREALF(x) (x.r)
    //#define CIMAGF(x) (x.i)
    #define FLOORF(x) (floorf(x))
  #else
    #include <algorithm>
    #include <complex>
    typedef std::complex<float> floatcomplex;
    #define I std::complex<float>{0,1}
    #define CEXPF(x) (std::exp(x))
    #define CREALF(x) (x.r)
    #define CIMAGF(x) (x.i)
    #define MAX(x,y) (std::max(x,y))
    #define MIN(x,y) (std::min(x,y))
    #define FLOORF(x) (std::floor(x))
  #endif
#endif

struct parameters {
  long Nx;
  long Ny;
  float dx;
  float dy;
  float dz;
  long iz_start;
  long iz_end;
  unsigned char xSymmetry;
  unsigned char ySymmetry;
  float taperPerStep;
  float twistPerStep;
  float d;
  float n_0;
  floatcomplex *n_in;
  long  Nx_n;
  long  Ny_n;
  long  Nz_n;
  float dz_n;
  floatcomplex *Efinal;
  floatcomplex *E1;
  floatcomplex *E2;
  floatcomplex *Eyx;
  floatcomplex *n_out;
  floatcomplex *b;
  float *multiplier;
  floatcomplex ax;
  floatcomplex ay;
  float rho_e;
  float RoC;
  float sinBendDirection;
  float cosBendDirection;
  double precisePower;
  float precisePowerDiff;
  float EfieldPower;
};


void solve_bpm(int io, double *tdft, FILE *fp) {
    // HDF5ファイルの作成
    // 関数から?(fp の入替え?)
    hid_t file_id;
    // local
    hid_t group_id, dataset_id, dataspace_id, memspace_id;
    herr_t status;

    double fmax[] = {0, 0};
    char str[BUFSIZ];
    int converged = 0;

	// 加工
    //ここから
    struct parameters P_var;
    struct parameters *P = &P_var;
    P->Nx = Nx;
    P->Ny = Ny;
    P->dx = 0.001; // x軸の区切りデータの1step幅
    P->dy = 0.001; // y軸の区切りデータの1step幅
    P->dz = 0.001; // z軸の区切りデータの1step幅
    P->iz_start = 0; // z軸の区切りデータ先頭(index)
    P->iz_end = Nz; // z軸の区切りデータ末尾(index)
    //enum symmetry : unsignate char
    //  enumeration
    //    NoSymmetry   (0)
    //    Symmetry     (1)
    //    AntiSymmetry (2)
    //  end
    //end
    P->xSymmetry = 0; //NULL; //???(片側設定?)(enum値[0-2])(Array[1,1]?)
    P->ySymmetry = 0; //NULL; //???(片側設定?)(enum値[0-2])(Array[1,1]?)
    P->d = 0.0; //???
    P->n_0 = 1.41; //[] reference refractive index
    P->n_in = Zin; //? //(floatcomplex *)mxGetData(mxGetField(prhs[1],0,"n_mat"));
    //mwSize nDims = 2
    //mwSize const *dimPtr = mxGetDimensions(mxGetField(prhs[1],0,"n_mat"));
    // 違いが分からない(処理適用範囲?)
    P->Nx_n = Nx; // 
    P->Ny_n = Ny; // 
    P->Nz_n = Nz; // nDims > 2? (long)dimPtr[2]: 1;
    P->dz_n = Nz; // ?
    //固定値?(ガウシアンビームの設定?)
    //P.Lz = 5e-3;
    //P.taperScaling = 0.15;
    //P.twistRate = 2*pi/P.Lz;
    P->taperPerStep = 0.0; //*(float *)mxGetData(mxGetField(prhs[1],0,"taperPerStep"));
    P->twistPerStep = 0.0; //*(float *)mxGetData(mxGetField(prhs[1],0,"twistPerStep"));
    P->rho_e = 0.0; // ?(Array)
    P->RoC = 0.0; // ?(Array)
    P->sinBendDirection = 0.0; // sin(*(float *)mxGetData(mxGetField(prhs[1],0,"bendDirection"))/180*PI); // 向き(rad)
    P->cosBendDirection = 0.0; // cos(*(float *)mxGetData(mxGetField(prhs[1],0,"bendDirection"))/180*PI); // 向き(rad)
    //初期電界入力?
	//BPM-MATLAB ではコールバック関数の登録で対応
    P->E1 = NULL; // Input E field(Array)
    //結果格納?
    P->Efinal = NULL; // Output E field(Array)
    P->n_out = NULL; //Output refractive index(Array)
    //ガウシアンビームの出力値?
    P->precisePower = 10.0; // (float)mxGetScalar(mxGetField(prhs[1],0,"inputPrecisePower"));
    #ifndef __NVCC__
    P->E2 = (floatcomplex *)((P->iz_end - P->iz_start)%2? P->Efinal: malloc(P->Nx*P->Ny*sizeof(floatcomplex)));
    #endif
    P->multiplier = NULL; // Array of multiplier values to apply to the E field after each step, due to the edge absorber outside the main simulation window
	//P->ax = floatcomplex{0.0,0.0}; //*(floatcomplex *)mxGetData(mxGetField(prhs[1],0,"ax"));
    //P->ay = floatcomplex{0.0,0.0}; //*(floatcomplex *)mxGetData(mxGetField(prhs[1],0,"ay"));

    P->EfieldPower = 0;
    P->precisePowerDiff = 0;
    #ifdef __NVCC__
    int temp, nBlocks; gpuErrchk(cudaOccupancyMaxPotentialBlockSize(&nBlocks,&temp,&substep1a,0,0));
    dim3 blockDims(TILE_DIM,TILE_DIM,1);

    struct parameters *P_dev;
    struct debug D_var = {{0.0,0.0,0.0},{0,0,0}};
    struct debug *D = &D_var;
    struct debug *D_dev;
    createDeviceStructs(P,&P_dev,D,&D_dev);
    #else
    #ifdef _OPENMP
    bool useAllCPUs = false; //mxIsLogicalScalarTrue(mxGetField(prhs[1],0,"useAllCPUs"));
    long numThreads = useAllCPUs || omp_get_num_procs() == 1? omp_get_num_procs(): omp_get_num_procs()-1;
    #else
    long numThreads = 1;
    #endif
    P->b = (floatcomplex *)malloc(numThreads*MAX(P->Nx,P->Ny)*sizeof(floatcomplex));
    #ifdef _OPENMP
    #pragma omp parallel num_threads(useAllCPUs || omp_get_num_procs() == 1? omp_get_num_procs(): omp_get_num_procs()-1)
    #endif
    #endif

    // initial field
    initfield();

    // 温度配列の初期化
    //int Nx = 100, Ny = 100, Nz = 100;
    double alpha = 0.01;  // 熱拡散係数
    //double *T = (double *)malloc(Nx * Ny * Nz * sizeof(double));
    //double *P_loss = (double *)malloc(Nx * Ny * Nz * sizeof(double));
    //memset(T, 0, Nx * Ny * Nz * sizeof(double));
    //memset(P_loss, 0, Nx * Ny * Nz * sizeof(double));
    //NN
    double *T = (double *)malloc(NFreq2 * NN * sizeof(double));
    double *P_losses = (double *)malloc(NFreq2 * NN * sizeof(double));
    memset(T, 0, NFreq2 * NN * sizeof(double));
    memset(P_losses, 0, NFreq2 * NN * sizeof(double));

    // セルの幅（空間ステップ）を計算
    double Dx = Xn[Nx] - Xn[0] / Nx;
    double Dy = Yn[Ny] - Yn[0] / Ny;
    double Dz = Zn[Nz] - Zn[0] / Nz;
    sprintf(str, "%.6f %.6f %.6f", Dx, Dy, Dz);
    fprintf(stdout, "%s\n", str);

    // HDF5ファイルの作成
    file_id = H5Fcreate(FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    // time step iteration
    int itime;
    double t = 0;
    //double sigma = 1e3;      // 導電率 [S/m]（適宜変更）
    int material_id = 0; // 使用したい材料のIDを設定
    //double sigma = get_conductivity(material_id);
    //double epsr = get_relative_permittivity(material_id);
    //double amur = get_relative_permeability(material_id);
    double mu_double_prime = 1e-3; // 磁気損失係数 [H/m]（適宜変更）
    //double frequency = 200e12;      // 周波数 [Hz]（適宜変更）
    //double omega = 2 * M_PI * frequency; // 角周波数 [rad/s]
    //for (itime = 0; itime <= Solver.maxiter; itime++) {
    //Step 時間の対応ではなく、Z軸の変化を元に対応していく。
    for(long iz=P->iz_start; iz<P->iz_end; iz++) {
        // update H
        t += 0.5 * Dt;
        updateHx(t);
        updateHy(t);
        updateHz(t);

        // ABC H
        if      (iABC == 0) {
            murH(numMurHx, fMurHx, Hx);
            murH(numMurHy, fMurHy, Hy);
            murH(numMurHz, fMurHz, Hz);
        }
        else if (iABC == 1) {
            pmlHx();
            pmlHy();
            pmlHz();
        }

        // PBC H
        if (PBCx) {
            pbcx();
        }
        if (PBCy) {
            pbcy();
        }
        if (PBCz) {
            pbcz();
        }

        // update E
        t += 0.5 * Dt;
        updateEx(t);
        updateEy(t);
        updateEz(t);

        // dispersion E
        if (numDispersionEx) {
            dispersionEx(t);
        }
        if (numDispersionEy) {
            dispersionEy(t);
        }
        if (numDispersionEz) {
            dispersionEz(t);
        }

        // ABC E
        if      (iABC == 1) {
            pmlEx();
            pmlEy();
            pmlEz();
        }

        // feed
        if (NFeed) {
            efeed(itime);
        }

        // inductor
        if (NInductor) {
            eload();
        }

        // point
        if (NPoint) {
            vpoint(itime);
        }

        // DFT
        //const double t0 = cputime();
        //dftNear3d(itime);
        //*tdft += cputime() - t0;

        #ifdef __NVCC__
            substep1a<<<nBlocks, blockDims>>>(P_dev); // xy -> yx
            substep1b<<<nBlocks, blockDims>>>(P_dev); // yx -> yx
            substep2a<<<nBlocks, blockDims>>>(P_dev); // yx -> xy
            substep2b<<<nBlocks, blockDims>>>(P_dev); // xy -> xy
            applyMultiplier<<<nBlocks, blockDims>>>(P_dev,iz,D_dev); // xy -> xy
        #else
            substep1a(P);
            substep1b(P);
            substep2a(P);
            substep2b(P);
            applyMultiplier(P,iz,NULL);
        #endif

        #ifdef _OPENMP
        #pragma omp master
        #endif
        {
            #ifdef __NVCC__
            if(iz+1 < P->iz_end) swapEPointers<<<1,1>>>(P_dev,iz);
            updatePrecisePower<<<1,1>>>(P_dev);
            gpuErrchk(cudaDeviceSynchronize()); // Wait until all kernels have finished
            #else
            if(iz+1 < P->iz_end) swapEPointers(P,iz);
            updatePrecisePower(P);
            #endif
        }
        #ifdef _OPENMP
        #pragma omp barrier
        #endif

        // average and convergence
        if ((itime % Solver.nout == 0) || (itime == Solver.maxiter)) {
            // average
            double fsum[2];
            average(fsum);

            // average (post)
            Eiter[Niter] = fsum[0];
            Hiter[Niter] = fsum[1];
            Niter++;

            // monitor
            if (io) {
                sprintf(str, "%7d %.6f %.6f", itime, fsum[0], fsum[1]);
                fprintf(fp,     "%s\n", str);
                fprintf(stdout, "%s\n", str);
                fflush(fp);
                fflush(stdout);

                // 各時間ステップごとにグループを作成
                char group_name[32];
                snprintf(group_name, sizeof(group_name), "/data%06d", itime);
                group_id = H5Gcreate(file_id, group_name, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

                // Eフィールドデータセットの作成と書き込み
                hsize_t e_dims[4] = {1, NFreq2, NN, 6};
                dataspace_id = H5Screate_simple(4, e_dims, NULL);
                dataset_id = H5Dcreate(group_id, "E", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

                // 書き込み用のメモリスペースを修正
                hsize_t mem_dims[1] = {6};
                memspace_id = H5Screate_simple(1, mem_dims, NULL);

                for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
                    int64_t n0 = ifreq * NN;
                    for (int nn = 0; nn < NN; nn++) {
                        double e_value[6] = {
                            cEx_r[n0 + nn], cEy_r[n0 + nn], cEz_r[n0 + nn],
                            cEx_i[n0 + nn], cEy_i[n0 + nn], cEz_i[n0 + nn]
                        };

                        hsize_t e_offset[4] = {0, ifreq, nn, 0};
                        hsize_t e_count[4] = {1, 1, 1, 6};
                        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, e_offset, NULL, e_count, NULL);

                        // 書き込み
                        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, e_value);
                        if (status < 0) {
                            fprintf(stderr, "Error writing E data at itime=%d, ifreq=%d, nn=%d\n", itime, ifreq, nn);
                        }
                    }
                }
                H5Dclose(dataset_id);
                H5Sclose(dataspace_id);

                // Hフィールドデータセットの作成と書き込み
                hsize_t h_dims[4] = {1, NFreq2, NN, 6};
                dataspace_id = H5Screate_simple(4, h_dims, NULL);
                dataset_id = H5Dcreate(group_id, "H", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

                for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
                    int64_t n0 = ifreq * NN;
                    for (int nn = 0; nn < NN; nn++) {
                        double h_value[6] = {
                            cHx_r[n0 + nn], cHy_r[n0 + nn], cHz_r[n0 + nn],
                            cHx_i[n0 + nn], cHy_i[n0 + nn], cHz_i[n0 + nn]
                        };

                        hsize_t h_offset[4] = {0, ifreq, nn, 0};
                        hsize_t h_count[4] = {1, 1, 1, 6};
                        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, h_offset, NULL, h_count, NULL);
                        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, h_value);
                        if (status < 0) {
                            fprintf(stderr, "Error writing H data at itime=%d, ifreq=%d, nn=%d\n", itime, ifreq, nn);
                        }
                    }
               }
                H5Dclose(dataset_id);
                H5Sclose(dataspace_id);

                // 複素数用のHDF5データ型を定義
                hid_t complex_datatype = H5Tcreate(H5T_COMPOUND, sizeof(d_complex_t));
                H5Tinsert(complex_datatype, "real", HOFFSET(d_complex_t, r), H5T_NATIVE_DOUBLE);
                H5Tinsert(complex_datatype, "imag", HOFFSET(d_complex_t, i), H5T_NATIVE_DOUBLE);

                // Surfaceフィールドデータセットの作成と書き込み
                hsize_t surf_dims[4] = {1, NFreq2, NN, 6};
                dataspace_id = H5Screate_simple(4, surf_dims, NULL);
                dataset_id = H5Dcreate(group_id, "Surface", complex_datatype, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

                for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
                    //int64_t surf0 = ifreq * NSurface;
                    for (int surf = 0; surf < NSurface; surf++) {
                        d_complex_t surf_value[6] = {
                            SurfaceEx[ifreq][surf], SurfaceEy[ifreq][surf], SurfaceEz[ifreq][surf],
                            SurfaceHx[ifreq][surf], SurfaceHy[ifreq][surf], SurfaceHz[ifreq][surf]
                        };

                        hsize_t surf_offset[4] = {0, ifreq, surf, 0};
                        hsize_t surf_count[4] = {1, 1, 1, 6};
                        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, surf_offset, NULL, surf_count, NULL);
                        status = H5Dwrite(dataset_id, complex_datatype, memspace_id, dataspace_id, H5P_DEFAULT, surf_value);
                        if (status < 0) {
                            fprintf(stderr, "Error writing H data at itime=%d, ifreq=%d, surf=%d\n", itime, ifreq, surf);
                        }
                    }
                }
                H5Dclose(dataset_id);
                H5Sclose(dataspace_id);

                // Pフィールドデータセットの作成と書き込み（仮の例）
                hsize_t p_dims[4] = {1, NFreq2, NN, 3};
                dataspace_id = H5Screate_simple(4, p_dims, NULL);
                dataset_id = H5Dcreate(group_id, "P", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

                // 書き込み用のメモリスペースを修正
                hsize_t mem_dims2[1] = {3};
                memspace_id = H5Screate_simple(1, mem_dims2, NULL);

                for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
                    int64_t n0 = ifreq * NN;
                    for (int nn = 0; nn < NN; nn++) {
                        double p_value[3] = {
                            cEx_r[n0 + nn] * cHy_r[n0 + nn] - cEy_r[n0 + nn] * cHx_r[n0 + nn],
                            cEy_r[n0 + nn] * cHz_r[n0 + nn] - cEz_r[n0 + nn] * cHy_r[n0 + nn],
                            cEz_r[n0 + nn] * cHx_r[n0 + nn] - cEx_r[n0 + nn] * cHz_r[n0 + nn]
                        };

                        hsize_t p_offset[4] = {0, ifreq, nn, 0};
                        hsize_t p_count[4] = {1, 1, 1, 3};
                        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, p_offset, NULL, p_count, NULL);
                        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, p_value);
                        if (status < 0) {
                            fprintf(stderr, "Error writing P data at itime=%d, ifreq=%d, nn=%d\n", itime, ifreq, nn);
                        }
                    }
                }
                H5Dclose(dataset_id);
                H5Sclose(dataspace_id);

                // 発熱量の計算
                hsize_t p_dims2[4] = {1, NFreq2, NN, 1};
                dataspace_id = H5Screate_simple(4, p_dims2, NULL);
                dataset_id = H5Dcreate(group_id, "P_loss", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

                // 書き込み用のメモリスペースを修正
                hsize_t mem_dims3[1] = {1};
                memspace_id = H5Screate_simple(1, mem_dims3, NULL);

                // 材料の特性設定
                for (int ifreq = 0; ifreq < NFreq2; ifreq++) {
                    int64_t n0 = ifreq * NN;
                    for (int nn = 0; nn < NN; nn++) {
                        // 発熱量密度の計算
                        double P_loss[1] = { P_losses[n0 + nn] };

                        hsize_t p_offset[4] = {0, ifreq, nn, 0};
                        hsize_t p_count[4] = {1, 1, 1, 1};
                        H5Sselect_hyperslab(dataspace_id, H5S_SELECT_SET, p_offset, NULL, p_count, NULL);
                        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, memspace_id, dataspace_id, H5P_DEFAULT, P_loss);
                        if (status < 0) {
                            fprintf(stderr, "Error writing P data at itime=%d, ifreq=%d, nn=%d\n", itime, ifreq, nn);
                        }
                    }
                }
                H5Dclose(dataset_id);
                H5Sclose(dataspace_id);

                // グループのクローズ
                H5Gclose(group_id);
            }

            // check convergence
            fmax[0] = MAX(fmax[0], fsum[0]);
            fmax[1] = MAX(fmax[1], fsum[1]);
            if ((fsum[0] < fmax[0] * Solver.converg) &&
                (fsum[1] < fmax[1] * Solver.converg)) {
                converged = 1;
                break;
            }
            
            // Niterを増加
            Niter++;
        }
    }
    // メモリの解放
    free(T);
    free(P_losses);

    // メモリスペース、データセットとデータスペースのクローズ
    status = H5Sclose(memspace_id);

    // result
    if (io) {
        sprintf(str, "    --- %s ---", (converged ? "converged" : "max steps"));
        fprintf(fp,     "%s\n", str);
        fprintf(stdout, "%s\n", str);
        fflush(fp);
        fflush(stdout);
    }

    // time steps
    Ntime = itime + converged;

    // メタデータの作成
    hid_t metadata_group_id = H5Gcreate(file_id, "/metadata", H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);

    // 時間に関するメタデータの書き込み(収束条件で終了時の対応)
    //int maxIter = MIN(Solver.maxiter, Ntime);
    //int maxNOut = MIN(Solver.nout, Niter);
    //double time_metadata[1] = {maxIter * Dt};
    double time_metadata[1] = {Solver.maxiter * Dt};
    //dataspace_id = H5Screate_simple(1, count, NULL);
    hsize_t time_count[1] = {1};
    dataspace_id = H5Screate_simple(1, time_count, NULL);
    dataset_id = H5Dcreate(metadata_group_id, "time", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, time_metadata);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

    // グリッドに関するメタデータの書き込み
    //double grid_metadata[3] = {Dx, Dy, Dz};
    //hsize_t grid_count[1] = {3};
    //dataspace_id = H5Screate_simple(1, grid_count, NULL);
    //dataset_id = H5Dcreate(metadata_group_id, "grid", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    //status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, grid_metadata);
    //H5Dclose(dataset_id);
    //H5Sclose(dataspace_id);

    // その他のメタデータの書き込み
/*
    // title, dt, source, fPlanewave, z0, Ni, Nj, Nk, N0, NN
    const char *title = Title;
    dataspace_id = H5Screate(H5S_SCALAR);
    dataset_id = H5Dcreate(metadata_group_id, "title", H5T_C_S1, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_C_S1, H5S_ALL, H5S_ALL, H5P_DEFAULT, title);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

    // 平面波が伝わるときのインピーダンスの値
    //double metadata_values[8] = {Dt, Planewave.z0, Ni, Nj, Nk, N0, NN};
    double metadata_values[8] = {Dt, 0.0, Ni, Nj, Nk, N0, NN};
    hsize_t metadata_count[1] = {8};
    dataspace_id = H5Screate_simple(1, metadata_count, NULL);
    dataset_id = H5Dcreate(metadata_group_id, "metadata_values", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, metadata_values);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

    // 配列に関するメタデータの書き込み (Xn, Yn, Zn, Freq1, Freq2)
    hsize_t array_count[1];
    double *arrays[] = {Xn, Yn, Zn, Freq1, Freq2};
    const char *array_names[] = {"Xn", "Yn", "Zn", "Freq1", "Freq2"};
    size_t array_sizes[] = {Nx + 1, Ny + 1, Nz + 1, NFreq1, NFreq2};

    for (int i = 0; i < 5; i++) {
        array_count[0] = array_sizes[i];
        dataspace_id = H5Screate_simple(1, array_count, NULL);
        dataset_id = H5Dcreate(metadata_group_id, array_names[i], H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, arrays[i]);
        H5Dclose(dataset_id);
        H5Sclose(dataspace_id);
    }
*/
    // Title
    hsize_t title_dims[1] = {256};
    dataspace_id = H5Screate_simple(1, title_dims, NULL);
    dataset_id = H5Dcreate(metadata_group_id, "Title", H5T_NATIVE_CHAR, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_CHAR, H5S_ALL, H5S_ALL, H5P_DEFAULT, Title);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

    // 各種整数型メタデータの書き込み
    struct {
        const char *name;
        void *value;
        hid_t type;
    } metadata[] = {
        {"Nx", &Nx, H5T_NATIVE_INT},
        {"Ny", &Ny, H5T_NATIVE_INT},
        {"Nz", &Nz, H5T_NATIVE_INT},
        {"Ni", &Ni, H5T_NATIVE_INT},
        {"Nj", &Nj, H5T_NATIVE_INT},
        {"Nk", &Nk, H5T_NATIVE_INT},
        {"N0", &N0, H5T_NATIVE_INT},
        {"NN", &NN, H5T_NATIVE_INT64},
        {"NFreq1", &NFreq1, H5T_NATIVE_INT},
        {"NFreq2", &NFreq2, H5T_NATIVE_INT},
        {"NFeed", &NFeed, H5T_NATIVE_INT},
        {"NPoint", &NPoint, H5T_NATIVE_INT},
        {"Niter", &Niter, H5T_NATIVE_INT},
        {"Ntime", &Ntime, H5T_NATIVE_INT},
        {"Solver_maxiter", &Solver.maxiter, H5T_NATIVE_INT},
        {"Solver_nout", &Solver.nout, H5T_NATIVE_INT},
        {"NGline", &NGline, H5T_NATIVE_INT},
        {"IPlanewave", &IPlanewave, H5T_NATIVE_INT}
    };

    for (int i = 0; i < sizeof(metadata) / sizeof(metadata[0]); i++) {
        dataspace_id = H5Screate(H5S_SCALAR);
        dataset_id = H5Dcreate(metadata_group_id, metadata[i].name, metadata[i].type, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, metadata[i].type, H5S_ALL, H5S_ALL, H5P_DEFAULT, metadata[i].value);
        H5Dclose(dataset_id);
        H5Sclose(dataspace_id);
    }

    // Dtの書き込み
    dataspace_id = H5Screate(H5S_SCALAR);
    dataset_id = H5Dcreate(metadata_group_id, "Dt", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Dt);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

    // Planewaveの書き込み
    dataspace_id = H5Screate(H5S_SCALAR);
    dataset_id = H5Dcreate(metadata_group_id, "Planewave", H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &Planewave);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

    // 配列データの書き込み
    struct {
        const char *name;
        double *data;
        size_t size;
    } arrays[] = {
        {"Xn", Xn, Nx + 1},
        {"Yn", Yn, Ny + 1},
        {"Zn", Zn, Nz + 1},
        {"Xc", Xc, Nx},
        {"Yc", Yc, Ny},
        {"Zc", Zc, Nz},
        {"Eiter", Eiter, Niter},
        {"Hiter", Hiter, Niter},
        {"VFeed", VFeed, NFeed * (Solver.maxiter + 1)},
        {"IFeed", IFeed, NFeed * (Solver.maxiter + 1)},
        {"VPoint", VPoint, NPoint * (Solver.maxiter + 1)},
        {"Freq1", Freq1, NFreq1},
        {"Freq2", Freq2, NFreq2},
        {"Gline", Gline, NGline * 2 * 3}
    };

    for (int i = 0; i < sizeof(arrays) / sizeof(arrays[0]); i++) {
        hsize_t array_dims[1] = {arrays[i].size};
        dataspace_id = H5Screate_simple(1, array_dims, NULL);
        dataset_id = H5Dcreate(metadata_group_id, arrays[i].name, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, arrays[i].data);
        H5Dclose(dataset_id);
        H5Sclose(dataspace_id);
    }
    
    // NSurfaceデータの書き込み
    dataspace_id = H5Screate(H5S_SCALAR);
    dataset_id = H5Dcreate(metadata_group_id, "NSurface", H5T_NATIVE_INT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, &NSurface);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);

    // Surfaceデータの書き込み
    // surface_t構造体に対応する複合データ型を定義
    hid_t memtype = H5Tcreate(H5T_COMPOUND, sizeof(surface_t));
    H5Tinsert(memtype, "nx", HOFFSET(surface_t, nx), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype, "ny", HOFFSET(surface_t, ny), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype, "nz", HOFFSET(surface_t, nz), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype, "x", HOFFSET(surface_t, x), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype, "y", HOFFSET(surface_t, y), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype, "z", HOFFSET(surface_t, z), H5T_NATIVE_DOUBLE);
    H5Tinsert(memtype, "ds", HOFFSET(surface_t, ds), H5T_NATIVE_DOUBLE);

    hsize_t surface_dims[1] = {NSurface};
    dataspace_id = H5Screate_simple(1, surface_dims, NULL);
    dataset_id = H5Dcreate(metadata_group_id, "Surface", memtype, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
    status = H5Dwrite(dataset_id, memtype, H5S_ALL, H5S_ALL, H5P_DEFAULT, Surface);
    H5Dclose(dataset_id);
    H5Sclose(dataspace_id);
    H5Tclose(memtype);

    // メタデータグループのクローズ
    H5Gclose(metadata_group_id);

    status = H5Fclose(file_id);

    // free
    memfree2();

    #ifdef __NVCC__
    gpuErrchk(cudaDeviceSynchronize()); // Wait until all kernels have finished
    retrieveAndFreeDeviceStructs(P,P_dev,D,D_dev);
    //   printf("\nDebug: %.18e %.18e %.18e %llu %llu %llu\n          ",D->dbls[0],D->dbls[1],D->dbls[2],D->ulls[0],D->ulls[1],D->ulls[2]);
    #else
    //if(P->E1 != mxGetData(prhs[0]) && P->E1 != P->Efinal) free(P->E1); // Part of the reason for checking this is to properly handle ctrl-c cases
    free(P->b);
    #endif
    double *outputPrecisePowerPtr = NULL; //(double *)mxGetData(plhs[2] = mxCreateDoubleMatrix(1,1,mxREAL));
    *outputPrecisePowerPtr = P->precisePower;
    return;
}


