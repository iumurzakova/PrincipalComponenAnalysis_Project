#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

// Function to calculate the mean of each feature
vector<double> calculateMean(const vector<vector<double>>& data) {
    vector<double> mean(data[0].size(), 0.0);
    int numSamples = data.size();
    
    

    for (int i = 0; i < numSamples; ++i) {
        for (int j = 0; j < data[i].size(); ++j) {
            mean[j] += data[i][j];
        }
    }
    for (int j = 0; j < mean.size(); ++j) {
        mean[j] /= numSamples;
    }
    return mean;
}

// Function to subtract the mean from each feature and return the centered data
vector<vector<double>> subtractMean(const vector<vector<double>>& data, const vector<double>& mean) {
    vector<vector<double>> centeredData = data;
    for (int i = 0; i < centeredData.size(); ++i) {
        for (int j = 0; j < centeredData[i].size(); ++j) {
            centeredData[i][j] -= mean[j];
        }
    }
    return centeredData;
}

// Function to compute the covariance matrix
vector<vector<double>> computeCovarianceMatrix(const vector<vector<double>>& data) {
    int numSamples = data.size();
    int numFeatures = data[0].size();
    
    

    // Calculate mean of each feature
    vector<double> mean = calculateMean(data);

    // Subtract the mean from each feature
    vector<vector<double>> centeredData = subtractMean(data, mean);

    // Compute the covariance matrix
    vector<vector<double>> covarianceMatrix(numFeatures, vector<double>(numFeatures, 0.0));
    for (int i = 0; i < numFeatures; ++i) {
        for (int j = 0; j < numFeatures; ++j) {
            double sum = 0.0;
            for (int k = 0; k < numSamples; ++k) {
                sum += centeredData[k][i] * centeredData[k][j];
            }
            covarianceMatrix[i][j] = sum / numSamples; // Unbiased estimate
        }
    }

    return covarianceMatrix;
}

// Function to find eigenvalues of a 2x2 matrix
void findEigenvalues(double a, double b, double c, double d, double& lambda1, double& lambda2) {
    double trace = a + d;
    double determinant = a * d - b * c;
    double discriminant = trace * trace - 4 * determinant;

    if (discriminant >= 0) { // Real eigenvalues
        lambda1 = (trace + std::sqrt(discriminant)) / 2.0;
        lambda2 = (trace - std::sqrt(discriminant)) / 2.0;
    } else {
        lambda1 = trace / 2.0;
        lambda2 = trace / 2.0; // For simplicity, we set both to the same in complex cases
    }
}

void findEigenvector(double a, double b, double c, double d, double lambda, double& vecX, double& vecY) {
    // Assuming matrix:
    // | a  b |
    // | c  d |
    // Eigenvector for lambda is (x, y) where:
    // (a - lambda)x + by = 0
    // cx + (d - lambda)y = 0
    if (a - lambda != 0) {
        vecY = 1;
        vecX =  - (b * vecY) / (a - lambda) ; // Solves for y
    } else if (b != 0) {
        vecX = 1;
        vecY = - (c * vecX) / (d - lambda) ; // Solves
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        for x
    } else {
        vecY = 1;
        vecX = 0;
    }
}
// Function to perform eigenvalue decomposition
pair<vector<double>, vector<vector<double>>> eigenDecomposition(const vector<vector<double>>& matrix) {
    int n = matrix.size();
    vector<vector<double>> eigenvectors(n, vector<double>(n, 0.0));
    vector<double> eigenvalues(n, 0.0);

    // Perform eigenvalue decomposition for a 2x2 matrix
    findEigenvalues(matrix[0][0], matrix[0][1], matrix[1][0], matrix[1][1], eigenvalues[0], eigenvalues[1]);

    // Find eigenvectors
    findEigenvector(matrix[0][0], matrix[0][1], matrix[1][0], matrix[1][1], eigenvalues[0], eigenvectors[0][0], eigenvectors[0][1]);
    findEigenvector(matrix[0][0], matrix[0][1], matrix[1][0], matrix[1][1], eigenvalues[1], eigenvectors[1][0], eigenvectors[1][1]);

    return make_pair(eigenvalues, eigenvectors);
}

    int main() {
    // Example usage
    vector<vector<double>> data = {{17.0,45.0},{18.0, 58.0},{18.0,59.0},{19.0,57.0},{19.0,65.0}};

    //Compute covariance matrix
    vector<vector<double>> covarianceMatrix = computeCovarianceMatrix(data);
// Output the covariance matrix
    cout << "Covariance Matrix:" << endl;
    for (const auto& row : covarianceMatrix) {
        for (double val : row) {
            cout << val << " ";
        }
        cout << endl;
    }

    // Perform eigenvalue decomposition
    auto result = eigenDecomposition(covarianceMatrix);
    vector<double> eigenvalues = result.first;
    vector<vector<double>> eigenvectors = result.second;

    // Output the eigenvalues
    cout << "\nEigenvalues:" << endl;
    for (double val : eigenvalues) {
        cout << val << " ";
    }
        cout << endl;

    // Output the eigenvectors
    cout << "\nEigenvectors:" << endl;
    for (const auto& row : eigenvectors) {
        for (double val : row) {
            cout << val << " ";
        }
        cout << endl;
    }

    // Calculate PCA results
    vector<vector<double>> pcaResults;
    for (const auto& sample : data) {
        vector<double> pcaResult;
        
        for (size_t i = 0; i < eigenvectors.size(); ++i) {
            double component = 0.0;
            
            for (size_t j = 0; j < eigenvectors[i].size(); ++j) {
                component += eigenvectors[i][j] * (sample[j] - calculateMean(data)[j]);
            }
             
        pcaResult.push_back(component);
        }
        
    pcaResults.push_back(pcaResult);
    }

    // Output the PCA results
    cout << "\nPCA Results:" << endl;
    for (const auto& row : pcaResults) {
        for (double val : row) {
            cout << val << " ";
        }
        cout << endl;
    }

    return 0;
    }
