#include <stdio.h>
#include <stdlib.h>

#define N 32
typedef float type;

int main(int argc, char const *argv[]) {

    FILE *fp;
    int i, j;
    char s;
    type **array;

    fp = fopen("testmatrix", "rb");

    array=(type **)malloc(N * sizeof(type *));
    *array=(type *)malloc(N * N * sizeof(type));
    for(i = 1; i < N; i++)    array[i] = array[i-1] + N;

    fread(&array[0][0], N * N * sizeof(type), 1, fp);  //注意不能是地址 array
    fclose(fp);

    fp = fopen("testmatrix.dat", "w");
    for(i = 0; i < N; i++){
        for(j = 0; j < N; j++) fprintf(fp,"% f ",array[i][j]);
        putc('\n',fp);
    }
    fclose(fp);
    printf("transform complete, press any key to continue...");
    s = getchar();
    return 0;
}
