#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <sys/time.h>

class LaplaceSerial {
private:
    std::vector<std::vector<double> > old_array;
    std::vector<std::vector<double> > new_array;
    int rows, cols;
    int max_its;
    double tolerance;

public:
    LaplaceSerial(int r, int c, int max_iter, double tol) {
        rows = r;
        cols = c;
        max_its = max_iter;
        tolerance = tol;
        
        old_array.resize(rows);
        new_array.resize(rows);
        for (int i = 0; i < rows; i++) {
            old_array[i].resize(cols, 0.0);
            new_array[i].resize(cols, 0.0);
        }
    }
    
    void inBoundaries() {
        
        for (int j = 0; j < cols; j++) {
            old_array[0][j] = 100.0;
            new_array[0][j] = 100.0;
        }
        
        for (int j = 0; j < cols; j++) {
            old_array[rows-1][j] = 0.0;
            new_array[rows-1][j] = 0.0;
        }
        
        for (int i = 0; i < rows; i++) {
            old_array[i][0] = 50.0;
            new_array[i][0] = 50.0;
        }
        
        for (int i = 0; i < rows; i++) {
            old_array[i][cols-1] = 75.0;
            new_array[i][cols-1] = 75.0;
        }
    }
    
    double solveIteration() {
        double max_diff = 0.0;
        
        for (int i = 1; i < rows-1; i++) {
            for (int j = 1; j < cols-1; j++) {
                new_array[i][j] = 0.25 * (old_array[i-1][j] + old_array[i+1][j] + 
                                         old_array[i][j-1] + old_array[i][j+1]);
                
                double diff = fabs(new_array[i][j] - old_array[i][j]);
                if (diff > max_diff) {
                    max_diff = diff;
                }
            }
        }
        
        return max_diff;
    }
    
    void copyArrays() {
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                old_array[i][j] = new_array[i][j];
            }
        }
    }
    
    void printArray() {
        if (rows <= 10 && cols <= 10) {
            std::cout << "Elements in the Array:" << std::endl;
            for (int i = 0; i < rows; i++) {
                for (int j = 0; j < cols; j++) {
                    std::cout << std::setw(8) << std::fixed << std::setprecision(2) << new_array[i][j] << " ";
                }
                std::cout << std::endl;
            }
            std::cout << std::endl;
        } else {
            std::cout << "Array too large to print (size: " << rows << "x" << cols << ")" << std::endl;
        }
    }
    
    double getTime() {
        struct timeval tv;
        gettimeofday(&tv, NULL);
        return tv.tv_sec + tv.tv_usec / 1000000.0;
    }
    
    void solve() {
        std::cout << "Meow ^_^ getting ready ^_^" << std::endl;
        std::cout << "Grid size: " << rows << " x " << cols << std::endl;
        std::cout << "Maximum iterations: " << max_its << std::endl;
        std::cout << "Tolerance: " << tolerance << std::endl << std::endl;
        
        inBoundaries();
        
        std::cout << "Initial array:" << std::endl;
        printArray();
        
        double start_time = getTime();
        
        int iteration;
        double max_diff = 0.0;
        
        for (iteration = 0; iteration < max_its; iteration++) {
            max_diff = solveIteration();
            
            if (iteration % 100 == 0) {
                std::cout << "Iteration " << iteration << ", Max difference: " << max_diff << std::endl;
            }
            
            if (max_diff < tolerance) {
                std::cout << "Converged after " << iteration << " iterations" << std::endl;
                break;
            }
            
            copyArrays();
        }
        
        double end_time = getTime();
        double duration = (end_time - start_time) * 1000;
        
        std::cout << "Final array:" << std::endl;
        printArray();
        
        std::cout << "Serial execution time: " << duration << " milliseconds" << std::endl;
        std::cout << "Final max difference: " << max_diff << std::endl;
        std::cout << "Total iterations: " << iteration << std::endl;
        std::cout << "All done yeyy" << std::endl;
    }
};

int main() {
    int rows, cols, max_its;
    double tolerance;
    
    std::cout << "Hello to my SuperComputer!!!" << std::endl;
    std::cout << "Enter number of rows: ";
    std::cin >> rows;
    std::cout << "Enter number of columns: ";
    std::cin >> cols;
    std::cout << "Enter maximum iterations: ";
    std::cin >> max_its;
    std::cout << "Enter tolerance: ";
    std::cin >> tolerance;
    
    std::cout << std::endl;
    
    LaplaceSerial solver(rows, cols, max_its, tolerance);
    solver.solve();
    
    return 0;
}
