//
//  Sparse.swift
//  MyMatrix
//
//  Created by 酒井裕司 on 2025/03/07.
//
@testable import MyMatrix
import Accelerate
import Testing
// QR分解の実行と検証

struct Test3 {
    @Test func testSparseQRDecomposition() async throws{
        let rowIndices: [Int32] =    [ 0,  1, 1,  2]
        let columnIndices: [Int32] = [ 2,  0, 2,  1]
        let aValues: [Double] =      [10, 20, 5, 50]
        let yValues:[Double] = [30,35,100]
        let y2Values:[Double] = [300,350,1000]
        let expect:[Double] = [1,2,3]
        let expect2:[Double] = [10,20,30]
        let A = SparseConvertFromCoordinate(3, 3,
                                            4, 1,
                                            SparseAttributes_t(),
                                            rowIndices, columnIndices,
                                            aValues)
        defer {
            SparseCleanup(A)
        }
        print("epsilon = \(epsilon)")
        let x = solveSparse(A,yValues)
        print("answer is \(x)")
        #expect(almostEqualVectorsAsArray(x,expect))
        let org = sparseMulVec(A, x)
        #expect(almostEqualVectorsAsArray(org,yValues))
        let x2 = solveSparse(A,y2Values)
        print("answer is \(x2)")
        #expect (almostEqualVectorsAsArray(x2,expect2))
        let org2 = sparseMulVec(A, x2)
        #expect(almostEqualVectorsAsArray(org2,y2Values))
    }
    @Test func testSparseMultiVec() async{
        let rowCount = Int32(4)
        let columnCount = Int32(4)
        let blockCount = 4
        let blockSize = UInt8(1)
        let rowIndices: [Int32] = [0, 3, 0, 3]
        let columnIndices: [Int32] = [0, 0, 3, 3]
        let data:[Double] = [1.0, 4.0, 13.0, 16.0]
        let A = SparseConvertFromCoordinate(rowCount, columnCount,
                                            blockCount, blockSize,
                                            SparseAttributes_t(),
                                            rowIndices, columnIndices,
                                            data)
        defer {
            SparseCleanup(A)
        }
        let xValues:[Double] = [10.0, -1.0, -1.0, 10.0]
        let expected:[Double] = [ 140.0, 0.0, 0.0,  200.0 ]
        
        let res = sparseMulVec(A, xValues)
        #expect(almostEqualVectorsAsArray(res,expected))
    }
}

func almostEqualT<T:BinaryFloatingPoint>(_ a: T, _ b: T, epsilon: T = 0.001) -> Bool {
    let diff = abs(a - b)
    let relativeError = diff / max(abs(a), abs(b), 1.0)
    return relativeError < epsilon
}
func almostEqualVectorsAsArray<T:BinaryFloatingPoint> (_ a:[T],_ b:[T], epsilon: T = 0.001)  -> Bool {
    
    guard a.count == b.count else {
        return false
    }
    let innerP = zip(a,b).reduce(T(0.0)){$0 + $1.0 * $1.1}
    let aNorm = sqrt(a.reduce(T(0.0)){$0 + $1 * $1})
    let bNorm = sqrt(b.reduce(T(0.0)){$0 + $1 * $1})
    let prodN = aNorm * bNorm
    if aNorm == 0 && bNorm == 0 {
        return true
    } else if prodN == 0 {
        return false
    }
    let res = innerP/prodN
    if (T(1.0) - res) < epsilon {
        return true
    }
    return false
}

func createRandomSparseMatrix(rows: Int, cols: Int, density: Double) -> SparseMatrix_Double {
    var values: [Double] = []
    var rowIndices: [Int32] = []
    var colIndices: [Int32] = []

    for j in 0..<cols {
        for i in 0..<rows {
            if Double.random(in: 0...1) < density {
                values.append(Double.random(in: -1...1))
                rowIndices.append(Int32(i))
                colIndices.append(Int32(j))
            }
        }
    }
    let A = SparseConvertFromCoordinate(Int32(rows), Int32(cols),
                                        values.count, 1,
                                        SparseAttributes_t(),
                                        rowIndices, colIndices,
                                        values)

    
    return A
}
func createRondomVector(_ n:Int) -> [Double] {
    var vec:[Double] = []
    for _ in 0..<n{
        vec.append(Double.random(in: 0...1))
    }
    return vec
}
func convToSparse(_ a:[[Double]]) -> SparseMatrix_Double {
    let row = a.count
    let col = a[0].count
    var colIndices:[Int32] = []
    var rowIndices:[Int32] = []
    var values:[Double] = []
    
    for i in 0..<row {
        for j in 0..<col {
            if a[i][j] != 0.0 {
                rowIndices.append(Int32(i))
                colIndices.append(Int32(j))
                values.append(a[i][j])
            }
        }
    }
    let A = SparseConvertFromCoordinate(Int32(row), Int32(col),
                                        values.count, 1,
                                        SparseAttributes_t(),
                                        rowIndices, colIndices,
                                        values)
    return A
}
