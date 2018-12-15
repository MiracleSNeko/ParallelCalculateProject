#include <stdio.h>
#include <stdlib.h>

#define N 1024
typedef float type;

int main(int argc, char const *argv[]) {

    FILE *fp;
    int i, j;
    char s;
    type **array;

    fp = fopen("vector", "rb");

    array=(type **)malloc(N * sizeof(type *));
    *array=(type *)malloc(N * 1 * sizeof(type));
    for(i = 1; i < N; i++)    array[i] = array[i-1] + 1;

    fread(&array[0][0], N * 1 * sizeof(type), 1, fp);  //注意不能是地址 array
    fclose(fp);

    fp = fopen("vector.dat", "w");
    for(i = 0; i < N; i++){
        for(j = 0; j < 1; j++) fprintf(fp,"% f ",array[i][j]);
        putc('\n',fp);
    }
    fclose(fp);
    printf("transform complete, press any key to continue...");
    s = getchar();
    return 0;
}
