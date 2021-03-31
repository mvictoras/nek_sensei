#include "../Bridge.h"

#define ELEMS 256
#define DIM 8

int main( int argc, char** argv ) {

    MPI_Init(NULL, NULL);
    MPI_Comm comm = MPI_COMM_WORLD;
    int elems = ELEMS;
    int x_dim = DIM;
    int y_dim = DIM;
    int z_dim = 1;
    int t_dim = 0;
    char casename[] = "test_ns"; 

    double mesh_x[ELEMS][1][DIM][DIM];
    double mesh_y[ELEMS][1][DIM][DIM];
    double mesh_z[ELEMS][1][DIM][DIM];
    double vel_x[ELEMS][1][DIM][DIM];
    double vel_y[ELEMS][1][DIM][DIM];
    double vel_z[ELEMS][1][DIM][DIM];
    double pr[ELEMS][1][DIM][DIM];
    for(int elem=0; elem<ELEMS; elem++) {
        for(int z=0; z<1; z++) {
            for(int y=0; y<DIM; y++) {
                for(int x=0; x<DIM; x++) {
                    mesh_x[elem][z][y][x] = 0.1*(x+8*(elem%16));
                    mesh_y[elem][z][y][x] = 0.1*(y+8*(elem/16));
                    mesh_z[elem][z][y][x] = 0.0;
                    vel_x[elem][z][y][x] = 0.0;
                    vel_y[elem][z][y][x] = 0.0;
                    vel_z[elem][z][y][x] = 0.0;
                    pr[elem][z][y][x] = x+y+elem;
                }
            }
        }
    }

    sensei_bridge_initialize(&comm, casename, &elems,
            &x_dim, &y_dim, &z_dim, mesh_x[0][0][0], mesh_y[0][0][0], mesh_z[0][0][0],
            &x_dim, &y_dim, &z_dim, mesh_x[0][0][0], mesh_y[0][0][0], mesh_z[0][0][0],
            vel_x[0][0][0], vel_y[0][0][0], vel_z[0][0][0], pr[0][0][0], NULL, &t_dim);

    double time = 0.0;
    for( int i = 0; i < 100; i++ ) {
        sensei_bridge_update( &i, &time );
        time += 0.1;
    }

    sensei_bridge_finalize();
    MPI_Finalize();

    return 0;
}
