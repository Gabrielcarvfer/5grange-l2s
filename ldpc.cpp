#include <armaMex.hpp>
#include <string>
#include <limits>
#include <unordered_map>
#include <tuple>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    if (nlhs > 1) mexErrMsgTxt("Only one output is acepted");
    
    if (nrhs != 2) mexErrMsgTxt("This function expects 2 arguments");
    
    // Get LDPC data
    unsigned int ldpc_Z          = static_cast<unsigned int>(mxGetScalar(mxGetField(prhs[1], 0, "Z")));
    unsigned int ldpc_numParBits = static_cast<unsigned int>(mxGetScalar(mxGetField(prhs[1], 0, "numParBits")));
    unsigned int ldpc_numTotBits = static_cast<unsigned int>(mxGetScalar(mxGetField(prhs[1], 0, "numTotBits")));
    unsigned int ldpc_numInfBits = static_cast<unsigned int>(mxGetScalar(mxGetField(prhs[1], 0, "numInfBits")));
    unsigned int ldpc_iterations = static_cast<unsigned int>(mxGetScalar(mxGetField(prhs[1], 0, "iterations")));
    
    std::string decType(mxArrayToString(mxGetField(prhs[1], 0, "decType")));
    
    arma::mat ldpc_H = armaGetPr(mxGetField(prhs[1], 0, "H"));
    
    // Get the LLR data
    arma::mat llrVec = armaGetPr(prhs[0]);
    unsigned int numBlocks = llrVec.n_rows;
    
    // Create the output decoded matrix
    plhs[0] = armaCreateMxMatrix(numBlocks, ldpc_numInfBits);
    arma::mat decVec = armaGetPr(plhs[0]);
    
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    
    unsigned int numEntries = arma::accu(ldpc_H > 0.5);
    
    double minVal = std::numeric_limits<double>::min();
    
    arma::mat Qv(numBlocks, 2*ldpc_Z + llrVec.n_cols);
    Qv.cols(0, 2*ldpc_Z - 1).fill(0);
    Qv.cols(2*ldpc_Z, 2*ldpc_Z + llrVec.n_cols - 1) = llrVec;
    
    std::unordered_map<unsigned int, std::unordered_map<unsigned int, arma::vec>> Rcv;
    
    // Find all neighbouring variable nodes of all nodes
    std::vector<arma::uvec> neighbouring;
    for (unsigned int checkIdx = 0; checkIdx < ldpc_numParBits; ++checkIdx) {
        neighbouring.push_back(arma::find(ldpc_H.row(checkIdx) == 1));
    }
    
    // Choose decoder
    if (decType == "SPA") {
        
        // Loop max number of algorithm iterations
        for (unsigned int ldpcCurIter = 0; ldpcCurIter < ldpc_iterations; ++ldpcCurIter) {
            
            // Loop over all check nodes
            for (unsigned int checkIdx = 0; checkIdx < ldpc_numParBits; ++checkIdx) {
                
                // Find all neighbouring variable nodes of current check node
                arma::uvec nbVarNodes = neighbouring[checkIdx];
                
                // Tmp update llr
                arma::mat tmpLlr(numBlocks, nbVarNodes.n_elem);
                for (unsigned int idx = 0; idx < nbVarNodes.n_elem; ++idx) {
                    if (Rcv[checkIdx].count(nbVarNodes(idx)) > 0) {
                        tmpLlr.col(idx) = Qv.col(nbVarNodes(idx)) - Rcv[checkIdx].at(nbVarNodes(idx));
                    } else {
                        tmpLlr.col(idx) = Qv.col(nbVarNodes(idx));
                    }
                }
                
                // Compute S = (Smag, Ssign)
                arma::vec Smag = arma::sum(-arma::log(minVal + arma::tanh(arma::abs(tmpLlr)/2)), 1);
                
                // Count number of negative elements
                arma::vec Ssign(numBlocks);
                arma::uvec tmp = arma::sum(tmpLlr < 0, 1);
                for (unsigned int idx = 0; idx < numBlocks; ++idx) {
                    Ssign(idx) = (tmp(idx) % 2)?-1:1;
                }
                
                // Loop all neighbouring variable nodes
                for (unsigned int varIter = 0; varIter < nbVarNodes.n_elem; ++varIter) {
                    unsigned int varIdx = nbVarNodes(varIter);
                    
                    arma::vec Qtmp(numBlocks);
                    if (Rcv[checkIdx].count(varIdx) > 0) {
                        Qtmp = Qv.col(varIdx) - Rcv[checkIdx].at(varIdx);
                    } else {
                        Qtmp = Qv.col(varIdx);
                    }
                    
                    arma::vec QtmpMag = -arma::log(minVal + arma::tanh(arma::abs(Qtmp)/2));
                    
                    // Note: +minVal in order to deal with llr = 0;
                    // implementation can be improved
                    arma::vec QtmpSign = sign(Qtmp + minVal);
                    
                    // Update message passing matrix
                    // From reference: Rcv = phi^-1(S-phi(Qtmp))
                    Rcv[checkIdx][varIdx] = Ssign % QtmpSign % (-arma::log(minVal + arma::tanh(arma::abs(Smag - QtmpMag)/2)));
                                       
                    // Update Qv. From reference: Qv = Qtmp + Rcv
                    Qv.col(varIdx) = Qtmp + Rcv[checkIdx].at(varIdx);
                }
            }
        }
        
        // Convert Qv to decoded bits
//         arma::mat decVec(numBlocks, ldpc_numInfBits, fill::zeros);
        decVec(arma::find(Qv.cols(0, ldpc_numInfBits - 1) < 0)).fill(1.0);
//         plhs[0] = armaCreateMxMatrix(decVec.n_rows, decVec.n_cols);
//         armaSetPr(plhs[0], decVec);
    }
    else {
        mexErrMsgTxt("Error: Unknown decoding algorithm LDPC.decType\n");
    }
}