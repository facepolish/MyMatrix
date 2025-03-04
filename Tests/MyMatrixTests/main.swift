//
//  main.swift
//  MyMatrix
//
//  Created by 酒井裕司 on 2025/02/23.
//
@testable import MyMatrix
import Testing
func almostEqual(_ a: Float, _ b: Float, epsilon: Float = 0.001) -> Bool {
    let diff = abs(a - b)
    let relativeError = diff / max(abs(a), abs(b), 1.0)
    return relativeError < epsilon
}
func almostEqualMatrix (_ a:Matrix,_ b:Matrix) -> Bool {
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
func almostEqualVectors (_ a:Vector,_ b:Vector, epsilon: Float = 0.001) throws -> Bool {
    
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
func isUpperTriangularMatrix(_ a: Matrix, epsilon: Float = 0.001) -> Bool {
    let mat: [[Float]] = a.get()
    for i in 0..<a.row {
        for j in 0..<i {
            if !almostEqual(mat[i][j], 0.0, epsilon: epsilon) {
                return false
            }
        }
    }
    return true
}
func isAlmostI(_ a: Matrix, epsilon: Float = 0.001) -> Bool {
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
func randMatrixRowCol (row:Int,col:Int) throws -> Matrix {
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
func randMatrix(_ rank:Int) throws -> Matrix {
    return try randMatrixRowCol(row:rank,col:rank)
}


struct Test {
    // Check QR decomposition
    func testQR(_ m:Matrix) throws {
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
 


 
