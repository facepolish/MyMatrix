//
//  Sparse.swift
//  MyMatrix
//
//  Created by 酒井裕司 on 2025/03/07.
//

import Accelerate
func solveSparse (_ m:SparseMatrix_Double,_ vec:[Double]) -> [Double] {
    let n = vec.count
    let factorization = sparseQRDecomposition(m)
    defer {
        SparseCleanup(factorization)
    }
    let m = vec.count
    var bValues = vec
    let xValues = [Double](unsafeUninitializedCapacity: vec.count) {
        xPtr, count in
        bValues.withUnsafeMutableBufferPointer { bPtr in
            let b = DenseVector_Double(
                count: Int32(n),
                data: bPtr.baseAddress!
            )
            let x = DenseVector_Double(
                count: Int32(m),
                data: xPtr.baseAddress!
            )
            SparseSolve(factorization, b, x)
        }
        count = vec.count
    }
    return xValues
}


func sparseQRDecomposition(_ a:SparseMatrix_Double) -> SparseOpaqueFactorization_Double{
    let symbolicOptions = SparseSymbolicFactorOptions(
        control: SparseDefaultControl,
        orderMethod: SparseOrderDefault,
        order: nil,
        ignoreRowsAndColumns: nil,
        malloc: { malloc($0) },
        free: { free($0) },
        reportError: nil)
    let numericOptions = SparseNumericFactorOptions()
    let factorization = SparseFactor(SparseFactorizationQR, a,
                                     symbolicOptions,
                                     numericOptions)
                                     
    return factorization
}

 func sparseMulVec (_ a:SparseMatrix_Double, _ b:[Double]) -> [Double] {
     let n = Int32(b.count)
     var xValues = b
     let yValues = [Double](unsafeUninitializedCapacity: xValues.count) {
         resultBuffer, count in
         xValues.withUnsafeMutableBufferPointer { denseMatrixPtr in
             let X = DenseVector_Double(count: n,
                                       data: denseMatrixPtr.baseAddress!)
             let Y = DenseVector_Double(count: n,
                                       data: resultBuffer.baseAddress!)
             
             SparseMultiply(a, X, Y)
         }
         count = xValues.count
     }
     return yValues
 }

