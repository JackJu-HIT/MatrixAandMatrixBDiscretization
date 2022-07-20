/**
 * @Function:Matrix A and Matrix B Discretization
 * @Date:2022-07-20 19:46:01
 * @Author:juchunyu
 * @Last modified:juchunyu
 */
#include <stdio.h>
#include <math.h>

#define SIZE  4  
#define T     0.02 
int main(){

       double MATRIX_A[SIZE][SIZE]                = {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16};
       double MATRIX_B[SIZE]                      = {1,2,3,4};
       double MATRIX_I[SIZE][SIZE]                = {1,0,0,0,0,1,0,0,0,0,1,0,0,0,0,1};
       double MATRIX_I_Multi_A[SIZE][SIZE]        = {0};
       double MATRIX_I_Sub_A[SIZE][SIZE]          = {0};
       double MATRIX_I_Sub_A_Inverse[SIZE][SIZE]  = {0};
       double MATRIX_A_Discretization[SIZE][SIZE] = {0};
       double MATRIX_B_Discretization[SIZE]       = {0};


       // B_D = T * B
       for(int i = 0;i < SIZE;i++){
           MATRIX_B_Discretization[i] = T * MATRIX_B[i];
       }

       /* A_D = (I+0.5*A*T)*(I-0.5*A*T)^(-1) */
       
       //MATRIX_I_Multi_A
       for(int i = 0;i < SIZE;i++){
           for(int j = 0;j < SIZE;j++){
               MATRIX_I_Multi_A[i][j] = MATRIX_I[i][j] + 0.5 * MATRIX_A[i][j] * T;
           } 
       }

       //MATRIX_I_Sub_A
        for(int i = 0;i < SIZE;i++){
           for(int j = 0;j < SIZE;j++){
               MATRIX_I_Sub_A[i][j] = MATRIX_I[i][j] - 0.5 * MATRIX_A[i][j] * T;
           } 
        }

        /** (I-0.5*A)^(-1) **/
        double Q[SIZE][SIZE]              = {0};
        double R[SIZE][SIZE]              = {0};
        double v[SIZE]                    = {0};
        double QT[SIZE][SIZE]             = {0};
        double RInverse[SIZE][SIZE]       = {0};
        double m_in[SIZE][SIZE]           = {0};
        double m_inverse_[SIZE][SIZE]     = {0};
        double m_inverse[SIZE][SIZE]      = {0}; 
        int    length                     = SIZE;  
    
        /**进行QR分解**/
        for(int k = 0;k < length;k++){

            /*R(1:k-1,k) = Q(:,1:k-1)’ * A(:,k)*/
            if(k >= 1){
                for(int i = 0;i < k;i++){
                    for(int j = 0;j < length;j++){
                        R[i][k] += Q[j][i]*MATRIX_I_Sub_A[j][k];
                    }
                }
            }

            /*v = A(:,k) - Q(:,1:k-1) * R(1:k-1,k)*/
            for(int j = 0;j<length;j++){
                if(k < 1){
                v[j] = MATRIX_I_Sub_A[j][0];
                } else {
                        v[j] = MATRIX_I_Sub_A[j][k];
                        for(int g = 0;g < k;g++){
                            v[j] -= R[g][k]*Q[j][g];
                        }
                }
            }

            /*R(k,k) = norm(v)*/
            for(int i = 0;i < length;i++){
                R[k][k] += v[i]*v[i];
            }

            R[k][k] = sqrt(R[k][k]);

            /*Q(:,k) = v / R(k,k)*/
            for(int i = 0;i < length;i++){
                Q[i][k] = v[i]/R[k][k];
            }  
        }
    
        //求解Q的转置
        for(int i = 0;i < length;i++){
        for(int j = 0;j<length;j++){
            QT[j][i] = Q[i][j];
        }
        }

        /**求解R的逆矩阵***/

        //转置
        for(int i = 0;i < length;i++){
            for(int j = 0;j < length;j++){
                m_in[j][i] = R[i][j];
            }
        }
        
    
        for(int j = 0;j < length;j++){
            //求解对角线
            m_inverse_[j][j] = 1/m_in[j][j];
            for(int i = j+1;i < length;i++){
                    double temp = 0;

                    for(int k = j;k <= i-1;k++){
                        temp += m_in[i][k]*m_inverse_[k][j];
                    }

                    m_inverse_[i][j] = -temp/m_in[i][i];
            }
        }
    
        //tranpose
        for(int i = 0;i < length;i++){
            for(int j = 0;j < length;j++){
                m_inverse[j][i] = m_inverse_[i][j];  //逆矩阵
            }
        }
    
        /*A = QR => A^(-1) = R^(-1)*Q^T*/

        for(int i = 0;i < length;i++){
            for(int j = 0;j < length;j++){
                for(int k = 0;k < length;k++){
                    MATRIX_I_Sub_A_Inverse[i][j] += m_inverse[i][k]*QT[k][j]; 
                }
            }
        }
    
        // A_D = (I+0.5*A*T)*(I-0.5*A*T)^(-1)
        for(int i = 0;i < SIZE;i++){
            for(int j = 0;j < SIZE;j++){
                for(int k = 0;k < SIZE;k++){
                    MATRIX_A_Discretization[i][j] += MATRIX_I_Multi_A[i][k] * MATRIX_I_Sub_A_Inverse[k][j];
                }
            }
        }
    
        //打印输出结果
        printf("MATRIX_A_Discretization：\n");
        for(int i = 0;i < SIZE;++i){
            for(int j = 0;j < SIZE;++j)
            {
                printf("%2.4f    ",MATRIX_A_Discretization[i][j]);
            }
            printf("\n");
        }

        printf("MATRIX_B_Discretization：\n");
        for(int i = 0;i < SIZE;++i){
            printf("%2.4f    ",MATRIX_B_Discretization[i]);
            printf("\n");
        }
}