//
//  main.swift
//  MyMatrix
//
//  Created by 酒井裕司 on 2025/02/23.
//

@testable import MyMatrix
import Testing
func almostEqual<T:BinaryFloatingPoint>(_ a: T, _ b: T, epsilon: T = 0.001) -> Bool {
    let diff = abs(a - b)
    let relativeError = diff / max(abs(a), abs(b), 1.0)
    return relativeError < epsilon
}
func almostEqualMatrix<T:BinaryFloatingPoint> (_ a:Matrix<T>,_ b:Matrix<T>) -> Bool {
    guard a.row == b.row && a.col == b.col else {
        return false
    }
    for item in zip(a.flat,b.flat){
        if !almostEqual(item.0,item.1) {
            return false
        }
    }
    return true
}
func almostEqualVectors<T:BinaryFloatingPoint> (_ a:Vector<T>,_ b:Vector<T>, epsilon: T = 0.001) throws -> Bool {
    
    guard a.size == b.size else {
        return false
    }
    let diff = try a - b
    let relativeError = diff.norm() / a.norm()
    if relativeError > epsilon {
        print("vector 1 = \(a.flat)")
        print("vector 2 = \(b.flat)")
        print ("relative error is \(relativeError)")
        return false
    }
    return true
}
func isUpperTriangularMatrix<T:BinaryFloatingPoint>(_ a: Matrix<T>, epsilon: T = 0.001) -> Bool {
    let mat: [[T]] = a.get()
    for i in 0..<a.row {
        for j in 0..<i {
            if !almostEqual(mat[i][j], 0.0, epsilon: epsilon) {
                return false
            }
        }
    }
    return true
}
func isAlmostI<T:BinaryFloatingPoint>(_ a: Matrix<T>, epsilon: T = 0.001) -> Bool {
    guard a.row == a.col else {
        return false
    }
    let mat = a.get()
    for i in 0..<a.row {
        for j in 0..<a.col {
            if ((i == j) && !almostEqual(mat[i][j], 1.0, epsilon: epsilon)) ||
               ((i != j) && !almostEqual(mat[i][j], 0.0, epsilon: epsilon)) {
                return false
            }
        }
    }
    return true
}
func randMatrixRowCol (row:Int,col:Int) throws -> Matrix<Float> {
    var b:[[Float]]=[]
    for _ in 0..<col{
        var a:[Float]=[]
        for _ in 0..<row{
            a.append(Float.random(in:0...1))
        }
        b.append(a)
    }
    return try Matrix(b)
}
func randMatrix(_ rank:Int) throws -> Matrix<Float> {
    return try randMatrixRowCol(row:rank,col:rank)
}


struct Test {
    // Check QR decomposition
    func testQR<T:BinaryFloatingPoint>(_ m:Matrix<T>) throws {
        let orth = try m.orthnormal()
        let II = try orth.transpose() * orth
        #expect( isAlmostI(II))
        let mm = try m.qrDecomposition()
        let mul = try mm.q.transpose() * mm.q
        #expect( isAlmostI( mul ) )
        #expect( isUpperTriangularMatrix(mm.r) )
    }
    @Test func testMatrixMultiply() async throws {
        let a:[[Float]] = [[1,2,3],[4,5,6]]
        let b:[[Float]] = [[10,11,12,13,14],[15,16,17,18,19],[20,21,22,23,24]]
        let c:[[Float]] = [[100,106,112,118,124],[235,250,265,280,295]]
        let d:[Float] = [4,5,6]
        let e:[Float] = [1,4]
        let a_mat = try Matrix(a)
        let b_mat = try Matrix(b)
        let c_mat = try Matrix(c)
        let res_mat = try a_mat * b_mat
        #expect( almostEqualMatrix(res_mat,c_mat))
        let d_vec = Vector(d)
        let e_vec = a_mat.vector(1,orientation: false)
        #expect( try almostEqualVectors(e_vec,d_vec))
        let f_vec = a_mat.vector(0,orientation: true)
        #expect( try almostEqualVectors(f_vec,Vector(e)))
        
        let la_ret = multiplyMatrix(a_mat, b_mat)
        #expect( almostEqualMatrix(la_ret,c_mat))
        print("epsilon = \(epsilon)")
    }
    @Test func testMatrixVector() async throws {
        let checkMatrix3:[[Float]] = [[0,1,1],[1,0,1],[1,1,0]]
        let vec1Flat:[Float] = [1.0,1,1]
        let vec2Flat:[Float] = [2.0,2,2]
        let vec2 = Vector(vec2Flat)
        let mat1 = try Matrix(checkMatrix3)
        let vec1 = Vector(vec1Flat)
        let res = try mat1 * vec1
        #expect( try almostEqualVectors(res,vec2) )
        print (res)
    }
    @Test func testMatrix() async throws {
        // Check QR decomposition
        let checkMatrix3:[[Float]] = [[0,1,1],[1,0,1],[1,1,0]]
        let checkMatrix4:[[Float]] = [[0,1,1,1],[1,0,1,1],[1,1,0,1]]
        let m = try Matrix(checkMatrix3)
        try testQR(m)
        let m2 = try Matrix(checkMatrix4)
        try testQR(m2)
        
        print(" test with col ,row <20 ")
        for _ in 0..<10{
            let n = Int.random(in: 2..<50)
            let m = Int.random(in: 2..<50)
            let mat = try randMatrixRowCol(row: n,col: m)
            try testQR(mat)
        }
    }
    
    @Test func testEigenFunc() async throws {
        let winningPoints:[Double] = [-26,31,-11,  -33,42,-15,  -25,23,-4]
        let ret = eigenvaluesAndEigenvectors(matrix:winningPoints,order:3)
        #expect( ret != nil )
        print(ret!.eigenvalues)
    }
    @Test func testGooleMethod() async throws {
        let winningPoints:[[Float]] = [[1,2,2],[3,1,1],[2,1,2]]
        let resultEigen:Float = 5.0
        let winMatrix = try Matrix(winningPoints)
        let gglp = try winMatrix.googleVector()
        #expect(gglp != nil)
        if let ggl = gglp {
            print("google vector is ")
            print( ggl.eigenVec.get())
            print("eigen value is " + String(ggl.eigen)   )
            #expect(abs(ggl.eigen - resultEigen) < 0.00001)
        }
        let rndMatrix = try randMatrix(10)
        let winMatRndp = try rndMatrix.googleVector()
        #expect(winMatRndp != nil)
        if let winMatRnd = winMatRndp {
            print("google vector is ")
            print( winMatRnd.eigenVec.flat)
            print("eigen value is " + String(winMatRnd.eigen)   )
        }
        let matrixWithImaginaryEigen:[[Float]] = [[3,1],[-5,-1]]
        let nonSquare:[[Float]] = [[1,3,1,1],[3,2,3,1],[3,3,5,5]]
//        let paraMatrix:[[Float]] = [[1,1,1],[1,1,1],[1,1,1]]
        let negativeMatrix = try Matrix(matrixWithImaginaryEigen)
        #expect(throws:(Errors.ParameterError)){
            let _ = try negativeMatrix.googleVector()
        }
        let nonSqMatrix = try Matrix(nonSquare)
        #expect(throws:(Errors.ParameterError)){
            let _ = try nonSqMatrix.googleVector()?.eigen
        }
    }
    @Test func testEigen() async throws {
   //        let winningPoints:[[Float]] = [[1,2,2],[3,1,1],[2,1,2]]
        let matrixWithImaginaryEigen:[[Float]] = [[3,1],[-5,-1]]
        let winningPoints:[[Float]] = [[-26,-33,-25],[31,42,23],[-11,-15,-4]]
        let winMatrix = try Matrix(winningPoints)
        let eigens = try winMatrix.realEigen()
        let count = eigens.eigenValues.count
        for i in 0..<count {
            let vec = eigens.eigenVectors[i]
            let eigenValue = eigens.eigenValues[i]
            let newVec = try winMatrix * vec
            let newVec2 = eigenValue * vec
            let boolvalue:Bool = try almostEqualVectors(newVec,newVec2)
            #expect( boolvalue )
            if !boolvalue {
                print ("eigen value is \(eigenValue)")
                print ("eigen vector is \(vec.flat)")
            }
        }
        let matIm = try Matrix(matrixWithImaginaryEigen)
        let imEigen = try matIm.realEigen()
        #expect(imEigen.eigenValues.count == 0)
    }
}

import Testing
//@testable import YourModuleName // Replace with the name of your module
func testAlmostEqual<T:BinaryFloatingPoint>(_ a: T, _ b: T, epsilon: T = 0.001) -> Bool {
    let diff = abs(a - b)
    let relativeError = diff / max(abs(a), abs(b), 1.0)
    return relativeError < epsilon
}
func testAlmostEqualMatrix<T:BinaryFloatingPoint> (_ a:[T],_ b:[T]) -> Bool {
    guard a.count == b.count else {
        return false
    }
    for item in zip(a,b){
        if !testAlmostEqual(item.0,item.1) {
            return false
        }
    }
    return true
}

struct  Test2 {
    @Test func testDGEMMCase1() async {
        // 入力行列の定義
        let a: [Double] = [
            1.0, 2.0, 3.0, 4.0,
            5.0, 6.0, 7.0, 8.0,
            9.0, 10.0, 11.0, 12.0
        ]
        let b: [Double] = [
            1.0, 2.0, 3.0,
            4.0, 5.0, 6.0,
            7.0, 8.0, 9.0,
            10.0, 11.0, 12.0
        ]

        // 期待される計算結果

        let expected: [Double] = [
            70.0, 80.0, 90.0,
            158.0, 184.0, 210.0,
            246.0, 288.0, 330.0
        ]
        // 行列積の実施
        let result = DGEMM(a: a, b: b, m: 3, n: 3, k: 4)
        #expect(testAlmostEqualMatrix(result, expected))
    }
    @Test func testDGEMMCase2() async {

        // 異なるテストセットの定義

        let a: [Double] = [
            2.0, 0.0, 1.0, 3.0,
            4.0, 1.0, 0.0, 2.0,
            5.0, 2.0, 1.0, 0.0
        ]
        let b: [Double] = [
            1.0, 2.0, 1.0,
            0.0, 1.0, 4.0,
            3.0, 4.0, 0.0,
            2.0, 0.0, 1.0
        ]
        // 期待される結果

        let expected: [Double] = [
            11.0, 8.0, 5.0,
            8.0, 9.0, 10.0,
            8.0, 16.0, 13.0
        ]
        // 行列積の実施
        let result = DGEMM(a: a, b: b, m: 3, n: 3, k: 4)
        #expect(testAlmostEqualMatrix(result, expected))
    }

}
struct Test4 {
    @Test func invertMatrixTest() async {
        // メイン処理
        let B: [Double] = [
            1, 3,
            2, 4
        ]
        let II: [[Double]] = [
            [1, 0],
            [0, 1]
        ]
        let invA = invert(matrix: B)
        guard let MatA = colMajorToMatrix(B, row: 2, col: 2), // 変換が成功した場合
              let matInv = colMajorToMatrix(invA, row: 2, col: 2) else {
            return
        }
        // 行列の乗算を行う
        let ret = matrixMultiply(matInv, MatA)
        print(ret)
        
        // 行列が等しいかを比較
        #expect (compareMatrices(ret, II) )
    }
    let matrix: [[Double]] = [
        [1.0, 2.0, 3.0],
        [4.0, 5.0, 6.0],
        [7.0, 8.0, 9.0]
    ]

    let matrix2: [[Double]] = [
        [2.0, 2.0, 2.0, 2.0],
        [1.0, -1.0, 1.0, -1.0],
        [-1, 1, -1, 1]
    ]
    @Test func testSVDcases() async throws{
        try testSVD(matrix)
        try testSVD(matrix2)
        func testSVD(_ a: [[Double]]) throws {
            if let svdResult = SVD(matrix: a) {
                print("左特異ベクトル行列U:")
                svdResult.U.forEach { print($0) }
                print("\n特異値:")
                print(svdResult.S)
                print("\n右特異ベクトル行列VT:")
                svdResult.VT.forEach { print($0) }
                
                let sMat = try #require(convSvaluetoMatrix(svdResult.S, row: a.count, col: a[0].count),"S matrix should not be nil")
                let m1 = matrixMultiply(svdResult.U, sMat)
                let m2 = matrixMultiply(m1, svdResult.VT)
                #expect(compareMatrices(m2, a) )
                if isOrthNormMatrix(svdResult.VT) {
                    print("right matrix is OK")
                }
                if isOrthNormMatrix(svdResult.U) {
                    print("left matrix is OK")
                }

                }
        }
    }
}
func isOrthNormMatrix<T:BinaryFloatingPoint>(_ a:[[T]],epsilon:T = 0.0001) -> Bool {
    let row = a.count
    let col = a[0].count
    guard row == col else {return false}
    let n = row
    for i in 0..<n {
        for j in 0..<n {
            var val:T = 0.0
            for k in 0..<n{
                val += a[i][k] * a[k][j]
            }
            if i == j {
                if abs(1.0 - val) > epsilon {
                    return false
                }
            } else {
                if abs(val) > epsilon{
                    return false
                }
            }
        }
    }
    return true
}

// 行列の乗算を行う関数
func matrixMultiply(_ a: [[Double]], _ b: [[Double]]) -> [[Double]] {
    let row = a.count
    let col = b[0].count // bの列数を取得
    let kCount = a[0].count
    var result = [[Double]](repeating: [Double](repeating: 0.0, count: col), count: row)
    
    for i in 0..<row {
        for j in 0..<col {
            for k in 0..<kCount {
                result[i][j] += a[i][k] * b[k][j] // 行列の積の計算
            }
        }
    }
    return result
}

// 行列が等しいかを比較する関数
func compareMatrices<T: BinaryFloatingPoint>(_ a: [[T]], _ b: [[T]], epsilon: T = 0.0001) -> Bool {
    guard a.count == b.count, a[0].count == b[0].count else {
        return false
    }
    for i in 0..<a.count {
        for j in 0..<a[0].count {
            let diff = abs(a[i][j] - b[i][j])
            if diff > epsilon {
                return false
            }
        }
    }
    return true
}
struct testVector {
    @Test func test1() async throws{
        let vector1:[Double] = [1.0, 2.0, 3.0]
        let vector2:[Double] = [2.0, 2.0, 2.0]
        let vector3:[Double] = [3.0,3.0,3.0,3.0]
        let vec1 = Vector(vector1)
        let vec2 = Vector(vector2)
        let vec3 = Vector(vector3)
        let expectedSize:Double = 12
        
        let result = try vec1 * vec2
        #expect(testAlmostEqual(result,expectedSize))
        print (result)
        #expect(throws:Errors.VectorSizeMismatch){
            _ = try vec1 * vec3
        }
    }
}
