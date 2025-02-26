//
//  main.swift
//  MyMatrix
//
//  Created by 酒井裕司 on 2025/02/23.
//
@testable import MyMatrix
import Testing

func almostEqual (_ a:Float,_ b:Float) -> Bool{
    let epsilon:Float = 0.001
    
    let diff = abs( a - b)
    //    let diffdiff = epsilon - diff
    if diff < epsilon {
//        print("error is " + String(diff))
        return true
    }
    return false
}
func almostEqualVectors (_ a:Vector,_ b:Vector) throws -> Bool {
    let epsilon:Float = 0.1
    
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
func isUpperTriangularMatrix(_ a:Matrix) -> Bool {
    let mat:[[Float]] = a.get()
    for i in 0..<a.row {
        for j in 0..<i {
            if !almostEqual(mat[i][j],0.0){
                return false
            }
        }
    }
    return true
}
func isAlmostI (_ a:Matrix) -> Bool{
    guard a.row == a.col else{
        return false
    }
    let mat = a.get()
    for i in 0..<a.row {
        for j in 0..<a.col {
            if  (( i == j) && !almostEqual(mat[i][j],1.0)) || ((i != j) && !almostEqual(mat[i][j], 0.0)){
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
    @Test func testGoogleVector() async throws {
        let winningPoints:[[Float]] = [[1,2,2],[3,1,1],[2,1,2]]
//        let resultVec:Vector = Vector([1,1,1])
        let resultEigen:Float = 5.0
        let winMatrix = try Matrix(winningPoints)
        let ggl = try winMatrix.googleVector()
        print("google vector is ")
        print( ggl.eigenVec.flat)
        print("eigen value is " + String(ggl.eigen)   )
        #expect(abs(ggl.eigen - resultEigen) < 0.001)
        let rndMatrix = try randMatrix(10)
        let winMatRnd = try rndMatrix.googleVector()
        print("google vector is ")
        print( winMatRnd.eigenVec.flat)
        print("eigen value is " + String(winMatRnd.eigen)   )
        
    }
    
    @Test func testEigenFunc() async throws {
        let winningPoints:[Double] = [-26,31,-11,  -33,42,-15,  -25,23,-4]
        let ret = eigenvaluesAndEigenvectors(matrix:winningPoints,order:3)
        #expect( ret != nil )
        print(ret!.eigenvalues)
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
 


 
