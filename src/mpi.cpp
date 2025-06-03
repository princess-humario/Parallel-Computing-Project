#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <mpi.h>
#include <sys/time.h>

class LaplaceMPI {
private:
    std::vector<std::vector<double> > local_old;
    std::vector<std::vector<double> > local_new;
    std::vector<double> top_ghost;
    std::vector<double> bottom_ghost;
    std::vector<double> send_top;
    std::vector<double> send_bottom;
    
    int global_rows, global_cols;
    int local_rows, local_cols;
    int max_its;
    double tolerance;
    int rank, size;
    int start_row, end_row;

public:
    LaplaceMPI(int g_rows, int g_cols, int max_iter, double tol) {
        global_rows = g_rows;
        global_cols = g_cols;
        max_its = max_iter;
        tolerance = tol;
        
        MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        MPI_Comm_size(MPI_COMM_WORLD, &size);
        
        // Calculate local domain size
        local_rows = global_rows / size;
        int remainder = global_rows % size;
        
        // Distribute remainder rows to first processes
        if (rank < remainder) {
            local_rows++;
            start_row = rank * local_rows;
        } else {
            start_row = rank * local_rows + remainder;
        }
        
        end_row = start_row + local_rows - 1;
        local_cols = global_cols;
        
        // Allocate local arrays with ghost rows
        int total_local_rows = local_rows;
        if (rank > 0) total_local_rows++; // top ghost row
        if (rank < size - 1) total_local_rows++; // bottom ghost row
        
        local_old.resize(total_local_rows);
        local_new.resize(total_local_rows);
        for (int i = 0; i < total_local_rows; i++) {
            local_old[i].resize(local_cols, 0.0);
            local_new[i].resize(local_cols, 0.0);
        }
        
        // Ghost row buffers
        top_ghost.resize(local_cols);
        bottom_ghost.resize(local_cols);
        send_top.resize(local_cols);
        send_bottom.resize(local_cols);
    }
    
    void initBoundaries() {
        int ghost_offset = (rank > 0) ? 1 : 0;
        
        // Top boundary (global row 0)
        if (start_row == 0) {
            for (int j = 0; j < local_cols; j++) {
                local_old[ghost_offset][j] = 100.0;
                local_new[ghost_offset][j] = 100.0;
            }
        }
        
        // Bottom boundary (global row global_rows-1)
        if (end_row == global_rows - 1) {
            int last_row = local_rows - 1 + ghost_offset;
            for (int j = 0; j < local_cols; j++) {
                local_old[last_row][j] = 0.0;
                local_new[last_row][j] = 0.0;
            }
        }
        
        // Left and right boundaries
        for (int i = 0; i < local_old.size(); i++) {
            local_old[i][0] = 50.0;
            local_new[i][0] = 50.0;
            local_old[i][local_cols-1] = 75.0;
            local_new[i][local_cols-1] = 75.0;
        }
    }
    
    void exchangeGhostRows() {
        int ghost_offset = (rank > 0) ? 1 : 0;
        
        // Prepare data to send
        if (rank > 0) {
            for (int j = 0; j < local_cols; j++) {
                send_top[j] = local_old[ghost_offset][j];  // First real row
            }
        }
        
        if (rank < size - 1) {
            int last_real_row = local_rows - 1 + ghost_offset;
            for (int j = 0; j < local_cols; j++) {
                send_bottom[j] = local_old[last_real_row][j];  // Last real row
            }
        }
        
        // Exchange with neighbors
        MPI_Request requests[4];
        int req_count = 0;
        
        // Send to top neighbor, receive from top neighbor
        if (rank > 0) {
            MPI_Isend(&send_top[0], local_cols, MPI_DOUBLE, rank-1, 0, MPI_COMM_WORLD, &requests[req_count++]);
            MPI_Irecv(&top_ghost[0], local_cols, MPI_DOUBLE, rank-1, 1, MPI_COMM_WORLD, &requests[req_count++]);
        }
        
        // Send to bottom neighbor, receive from bottom neighbor
        if (rank < size - 1) {
            MPI_Isend(&send_bottom[0], local_cols, MPI_DOUBLE, rank+1, 1, MPI_COMM_WORLD, &requests[req_count++]);
            MPI_Irecv(&bottom_ghost[0], local_cols, MPI_DOUBLE, rank+1, 0, MPI_COMM_WORLD, &requests[req_count++]);
        }
        
        // Wait for all communications to complete
        MPI_Waitall(req_count, requests, MPI_STATUSES_IGNORE);
        
        // Copy received ghost data
        if (rank > 0) {
            for (int j = 0; j < local_cols; j++) {
                local_old[0][j] = top_ghost[j];  // Top ghost row
            }
        }
        
        if (rank < size - 1) {
            int bottom_ghost_row = local_rows + ghost_offset;
            for (int j = 0; j < local_cols; j++) {
                local_old[bottom_ghost_row][j] = bottom_ghost[j];  // Bottom ghost row
            }
        }
    }
    
    double solveIteration() {
        double local_max_diff = 0.0;
        int ghost_offset = (rank > 0) ? 1 : 0;
        
        // Process internal points only
        int start_i = ghost_offset;
        int end_i = local_rows - 1 + ghost_offset;
        
        // Skip boundary rows for computation
        if (start_row == 0) start_i++;  // Skip global top boundary
        if (end_row == global_rows - 1) end_i--;  // Skip global bottom boundary
        
        for (int i = start_i; i <= end_i; i++) {
            for (int j = 1; j < local_cols - 1; j++) {
                local_new[i][j] = 0.25 * (local_old[i-1][j] + local_old[i+1][j] + 
                                         local_old[i][j-1] + local_old[i][j+1]);
                
                double diff = fabs(local_new[i][j] - local_old[i][j]);
                if (diff > local_max_diff) {
                    local_max_diff = diff;
                }
            }
        }
        
        return local_max_diff;
    }
    
    void copyArrays() {
        for (int i = 0; i < local_old.size(); i++) {
            for (int j = 0; j < local_cols; j++) {
                local_old[i][j] = local_new[i][j];
            }
        }
    }
    
    void printArray() {
        if (global_rows <= 10 && global_cols <= 10) {
            // Gather all data to rank 0 for printing
            if (rank == 0) {
                std::vector<std::vector<double> > global_array(global_rows);
                for (int i = 0; i < global_rows; i++) {
                    global_array[i].resize(global_cols);
                }
                
                // Copy rank 0's data
                int ghost_offset = 0;
                for (int i = 0; i < local_rows; i++) {
                    for (int j = 0; j < local_cols; j++) {
                        global_array[start_row + i][j] = local_new[i + ghost_offset][j];
                    }
                }
                
                // Receive data from other ranks
                for (int r = 1; r < size; r++) {
                    int other_local_rows, other_start_row;
                    MPI_Recv(&other_local_rows, 1, MPI_INT, r, 100, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    MPI_Recv(&other_start_row, 1, MPI_INT, r, 101, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    
                    for (int i = 0; i < other_local_rows; i++) {
                        MPI_Recv(&global_array[other_start_row + i][0], local_cols, MPI_DOUBLE, r, 102 + i, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
                    }
                }
                
                // Print the global array
                std::cout << "Elements in the Array:" << std::endl;
                for (int i = 0; i < global_rows; i++) {
                    for (int j = 0; j < global_cols; j++) {
                        std::cout << std::setw(8) << std::fixed << std::setprecision(2) << global_array[i][j] << " ";
                    }
                    std::cout << std::endl;
                }
                std::cout << std::endl;
            } else {
                // Send local data to rank 0
                MPI_Send(&local_rows, 1, MPI_INT, 0, 100, MPI_COMM_WORLD);
                MPI_Send(&start_row, 1, MPI_INT, 0, 101, MPI_COMM_WORLD);
                
                int ghost_offset = 1;  // Other ranks have top ghost row
                for (int i = 0; i < local_rows; i++) {
                    MPI_Send(&local_new[i + ghost_offset][0], local_cols, MPI_DOUBLE, 0, 102 + i, MPI_COMM_WORLD);
                }
            }
        } else {
            if (rank == 0) {
                std::cout << "Array too large to print (size: " << global_rows << "x" << global_cols << ")" << std::endl;
            }
        }
    }
    
    double getTime() {
        return MPI_Wtime();
    }
    
    void solve() {
        if (rank == 0) {
            std::cout << "Hello from MPI Land! We have " << size << " processes working together!" << std::endl;
            std::cout << "Grid size: " << global_rows << " x " << global_cols << std::endl;
            std::cout << "Maximum iterations: " << max_its << std::endl;
            std::cout << "Tolerance: " << tolerance << std::endl;
            std::cout << "Number of MPI processes: " << size << std::endl << std::endl;
        }
        
        initBoundaries();
        
        if (rank == 0) {
            std::cout << "Initial array:" << std::endl;
        }
        printArray();
        
        double start_time = getTime();
        
        int iteration;
        double local_max_diff, global_max_diff = 0.0;
        
        for (iteration = 0; iteration < max_its; iteration++) {
            // Exchange ghost rows before computation
            exchangeGhostRows();
            
            // Compute local iteration
            local_max_diff = solveIteration();
            
            // Find global maximum difference
            MPI_Allreduce(&local_max_diff, &global_max_diff, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
            
            if (rank == 0 && iteration % 100 == 0) {
                std::cout << "Iteration " << iteration << ", Max difference: " << global_max_diff << std::endl;
            }
            
            if (global_max_diff < tolerance) {
                if (rank == 0) {
                    std::cout << "Converged after " << iteration << " iterations!" << std::endl;
                }
                break;
            }
            
            copyArrays();
        }
        
        double end_time = getTime();
        double duration = end_time - start_time;
        
        // Find maximum time across all processes
        double max_time;
        MPI_Reduce(&duration, &max_time, 1, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
        
        if (rank == 0) {
            std::cout << "Final array:" << std::endl;
        }
        printArray();
        
        if (rank == 0) {
            std::cout << "MPI execution time: " << max_time << " seconds" << std::endl;
            std::cout << "Final max difference: " << global_max_diff << std::endl;
            std::cout << "Total iterations: " << iteration << std::endl;
            std::cout << "MPI processes used: " << size << std::endl;
        }
    }
};

int main(int argc, char* argv[]) {
    MPI_Init(&argc, &argv);
    
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    int rows, cols, max_its;
    double tolerance;
    
    if (rank == 0) {
        std::cout << "Welcome to MPI Supercomputer version!" << std::endl;
        std::cout << "Enter number of rows: ";
        std::cin >> rows;
        std::cout << "Enter number of columns: ";
        std::cin >> cols;
        std::cout << "Enter maximum iterations: ";
        std::cin >> max_its;
        std::cout << "Enter tolerance: ";
        std::cin >> tolerance;
        std::cout << std::endl;
    }
    
    // Broadcast input parameters to all processes
    MPI_Bcast(&rows, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&cols, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&max_its, 1, MPI_INT, 0, MPI_COMM_WORLD);
    MPI_Bcast(&tolerance, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    LaplaceMPI solver(rows, cols, max_its, tolerance);
    solver.solve();
    
    MPI_Finalize();
    return 0;
}
