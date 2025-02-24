//
//  Matrix.swift
//  Matrix
//
//  Created by 酒井裕司 on 2025/02/03.
//

import simd
import Foundation
import Accelerate
let epsilon:Float = 0.0001

enum Errors: Error {
    case MatrixSizeMismatch
    case VectorSizeMismatch
    case NoEigenValues
}
//typealias Scalar = Float
//typealias Scalar4x4 = float4x4

/*
func eigenvaluesAndEigenvectors(matrix: [Double], order: Int) -> (eigenvalues: [(real: Double, imag: Double)], eigenvectors: [[Double]])? {
    guard matrix.count == order * order else { return nil }

    var eigenvalueReal = [Double](repeating: 0.0, count: order)
    var eigenvalueImaginary = [Double](repeating: 0.0, count: order)
    var eigenvector = [Double](repeating: 0.0, count: order * order)
    var varMatrix = matrix
    var workspaceQuery = [Double](repeating: 0.0, count: 1)
    var status = Int32(0)

    // Workspace query
    var orderInt32 = Int32(order)
    var orderInt32_2 = Int32(order)
    var lda = Int32(order)
    var ldvr = Int32(order)
    var lwork = Int32(-1)
    var N = "N"//Int8("N")//"N"
    var V = "V"//Int8("V")//"V"

    dgeev_(&N, &V, &orderInt32, &varMatrix, &lda, &eigenvalueReal, &eigenvalueImaginary, nil, &orderInt32_2, &eigenvector, &ldvr, &workspaceQuery, &lwork, &status)

    guard status == 0 else { return nil }

    // Actual computation
    lwork = Int32(workspaceQuery[0])
    var workspace = [Double](repeating: 0.0, count: Int(lwork))

    dgeev_(&N,&V, &orderInt32, &varMatrix, &lda, &eigenvalueReal, &eigenvalueImaginary, nil, &orderInt32_2, &eigenvector, &ldvr, &workspace, &lwork, &status)

    guard status == 0 else { return nil }

    let eigenvalues = zip(eigenvalueReal, eigenvalueImaginary).map { (real: $0.0, imag: $0.1) }
    let eigenvectors = stride(from: 0, to: eigenvector.count, by: order).map { Array(eigenvector[$0..<$0+order]) }

    return (eigenvalues: eigenvalues, eigenvectors: eigenvectors)
}
 */
func eigenvaluesAndEigenvectors(matrix: [Double], order: Int) -> (eigenvalues: [(real: Double, imag: Double)], eigenvectors: [[Double]])? {
    guard matrix.count == order * order else { return nil }

    var eigenvalueReal = [Double](repeating: 0.0, count: order)
    var eigenvalueImaginary = [Double](repeating: 0.0, count: order)
    var eigenvector = [Double](repeating: 0.0, count: order * order)
    var varMatrix = matrix
    var workspaceQuery = [Double](repeating: 0.0, count: 1)
    var status = Int32(0)

    // Workspace query
    var orderInt32 = Int32(order)
    var orderInt32_2 = Int32(order)

    var lda = Int32(order)
    var ldvr = Int32(order)
    var lwork = Int32(-1)

    let N = "N".utf8CString
    let V = "V".utf8CString

    _ = N.withUnsafeBufferPointer { nPtr in
        V.withUnsafeBufferPointer { vPtr in
            dgeev_(UnsafeMutablePointer(mutating: nPtr.baseAddress),
                   UnsafeMutablePointer(mutating: vPtr.baseAddress),
                   &orderInt32, &varMatrix, &lda, &eigenvalueReal, &eigenvalueImaginary,
                   nil, &orderInt32_2, &eigenvector, &ldvr, &workspaceQuery, &lwork, &status)
        }
    }

    guard status == 0 else { return nil }

    // Actual computation
    lwork = Int32(workspaceQuery[0])
    var workspace = [Double](repeating: 0.0, count: Int(lwork))

    _ = N.withUnsafeBufferPointer { nPtr in
        V.withUnsafeBufferPointer { vPtr in
            dgeev_(UnsafeMutablePointer(mutating: nPtr.baseAddress),
                   UnsafeMutablePointer(mutating: vPtr.baseAddress),
                   &orderInt32, &varMatrix, &lda, &eigenvalueReal, &eigenvalueImaginary,
                   nil, &orderInt32_2, &eigenvector, &ldvr, &workspace, &lwork, &status)
        }
    }

    guard status == 0 else { return nil }

    let eigenvalues = zip(eigenvalueReal, eigenvalueImaginary).map { (real: $0.0, imag: $0.1) }
    let eigenvectors = stride(from: 0, to: eigenvector.count, by: order).map { Array(eigenvector[$0..<$0+order]) }

    return (eigenvalues: eigenvalues, eigenvectors: eigenvectors)
}

func convTo4x4(_ a:[[Float]]) throws -> float4x4 {
    let row = a.count
    let col = a[0].count
    if row != col || row != 4 {
        throw Errors.MatrixSizeMismatch
    }
    var ret = float4x4()
    for i in 0..<4 {
        ret[i] = SIMD4(a[i][0], a[i][1], a[i][2], a[i][3])
    }
    return ret
}
public class Matrix {
    public var flat: [Float] = []
    public var row: Int
    public var col: Int
    
    public init(_ a: [[Float]]) throws {
        self.row = a.count
        self.col = a[0].count
        for i in 0..<row {
            if a[i].count == col {
                self.flat = self.flat + a[i]
            } else {
                throw Errors.MatrixSizeMismatch
            }
        }
    }
    public func get(_ row:Int,_ col:Int) -> Float {
        let offset = self.col * row + col
        return self.flat[offset]
    }

    public init(_ a: [Float], row: Int, col: Int) {
        self.flat = a
        self.row = row
        self.col = col
    }
    public init(_ a:Vector,orientation:Bool){//if orientation is true -> column vector
        self.flat = a.flat
        if orientation {
            self.row = a.size
            self.col = 1
        }else{
            self.row = 1
            self.col = a.size
        }
    }

    public init(_ f:[Vector]){// assume vecgors as col vectors
        let cols = f.count
        let rows = f[0].size
        let flat_ini:[Float] = []
        let flat = f.reduce(flat_ini){$0+$1.flat}
        let retMat = Matrix(flat,row:cols,col:rows).transpose()
        self.flat = retMat.flat
        self.row = retMat.row
        self.col = retMat.col
    }
    public func orthnormal() throws -> Matrix {
        var vecs:[Vector] = []
        let count = min(self.col,self.row)
        for i in 0..<count {
            let vec = self.vector(i, orientation: true)
            var difVec = Vector(0.0,size:vec.size)
            difVec = try vecs.reduce(difVec){try $0 - (try vec*$1) * $1}
            let retVec = try vec + difVec
            vecs.append(retVec.normalize())
        }
        return Matrix(vecs) //.transpose()
    }
    public func qrDecomposition() throws -> (q:Matrix,r:Matrix) {
        let q = try self.orthnormal()
        let r = try q.transpose() * self
        return (q,r)
    }
    //extract row or column vector from matrix. if orientation is true -> column vec
    public func vector(_ index:Int,orientation:Bool) -> Vector {
        var flat:[Float] = []
        if orientation {// Column vector
            flat = stride(from:index,to:self.flat.count,by:self.col).map{i in self.flat[i]}
        }else{
            let offset = index*self.row
            let end = offset+self.col
            flat = Array(self.flat[offset..<end])
        }
        return Vector(flat)
    }
    public static func + (lhs: Matrix, rhs: Matrix) throws -> Matrix {
        guard lhs.row == rhs.row && lhs.col == rhs.col else{
            throw Errors.MatrixSizeMismatch
        }
        let plus = zip(lhs.flat, rhs.flat).map { $0.0 + $0.1 }
        return Matrix(plus, row: rhs.row, col: rhs.col)
    }
    private func isNonNegativeMatrix() -> Bool {
        guard self.col == self.row else {
            return false
        }
        let size = self.row * self.col
        for i in 0..<size {
            if self.flat[i] < 0 {
                return false
            }
        }
        return true
    }
    public func mainDiagonal() -> [Float] {
        var ret:[Float] = []
        let count = min(self.row,self.col)
        for i in 0..<count {
            ret.append(self.get(i,i))
        }
        return ret
    }
    public func googleVector() throws -> (eigen:Float,eigenVec:Vector) {// If the matrix is no-negative google vector exist
        guard self.isNonNegativeMatrix() else {
            throw Errors.MatrixSizeMismatch
        }
        let eigens = try self.eigen()
        var maxVal:Float = 0
        var maxVec:Vector = eigens.eigenVectors[0]
        for i in 0..<eigens.eigenValues.count {
            if eigens.eigenValues[i] > maxVal {
                maxVal = eigens.eigenValues[i]
                maxVec = eigens.eigenVectors[i]
            }
            
        }
        return (maxVal,maxVec)
    }
    public func eigen() throws -> (eigenValues:[Float], eigenVectors:[Vector]){
        guard self.row == self.col else {
            throw Errors.MatrixSizeMismatch
        }
        let count:Int = self.row * self.col
        var flat = [Double](repeating:0,count:count)
        for k in 0..<count{
            let i:Int = k%self.row
            let j:Int = k/self.col
            flat[k] = Double(self.get(i,j))
        }
        let ret = eigenvaluesAndEigenvectors(matrix:flat,order:self.row)
        guard let result = ret else {
            throw Errors.NoEigenValues
        }
        var retVectors:[Vector] = []
        var eigens:[Float] = []
        
        for item in zip(result.eigenvalues,result.eigenvectors){
            if item.0.imag != 0 {
                eigens.append(0)
            }else {
                eigens.append(Float(item.0.real))
            }
            let vec = Vector(item.1.map{Float($0)})
            retVectors.append(vec)
        }
        return (eigens,retVectors)
    }
    /*
    public func eigen() throws -> (eigenValues:[Float], eigenVectors:[Vector]){
        guard self.row == self.col else {
            throw Errors.MatrixSizeMismatch
        }
        let epsilon:Float = 0.000001
        let maxCount = 10000
        let res = try self.qrDecomposition()
        var q = res.q
        var r = res.r
        var a = try r * q
        var x = q
        var lastEigen = a.mainDiagonal()
        var diff:Float = 0
        var count = 0
        repeat {
            let res = try a.qrDecomposition()
            q = res.q
            r = res.r
            a = try r * q
            x = try x * q
            let topEigen = a.mainDiagonal()
            let diffVec = try Vector(lastEigen) - Vector(topEigen)
            diff = diffVec.norm()
//            print ("diff norm is \(diff)")
            lastEigen = topEigen
            count += 1
            if count > maxCount {
                break
            }
        }while(diff > epsilon)
        var retVectors:[Vector] = []
        for i in 0..<self.col {
            retVectors.append(x.vector(i,orientation: true))
        }
        return (lastEigen,retVectors)
    }
     */
    public static func - (lhs: Matrix, rhs: Matrix) throws -> Matrix {
        guard lhs.row == rhs.row && lhs.col == rhs.col else{
            throw Errors.MatrixSizeMismatch
        }
        let minus = zip(lhs.flat, rhs.flat).map { $0.0 - $0.1 }
        return Matrix(minus, row: rhs.row, col: rhs.col)
    }
    public static func * (lhs:Matrix, rhs:Vector) throws -> Vector { // Assume vector as col vector
        guard lhs.col == rhs.size else {
            throw Errors.MatrixSizeMismatch
        }
        let matt = lhs.get()
        var flat:[Float] = []
        for i in 0..<lhs.row {
            var val:Float = 0
            for j in 0..<lhs.col {
                val = val + matt[i][j] * rhs.flat[j]
            }
            flat.append(val)
        }
        return Vector(flat)
    }
    public static func * (lhs: Matrix, rhs: Matrix) throws -> Matrix {
        guard lhs.col == rhs.row else{
            throw Errors.MatrixSizeMismatch
        }
        var mul = [Float](repeating: 0.0, count: lhs.row * rhs.col)
        vDSP_mmul(lhs.flat, 1, rhs.flat, 1, &mul, 1, vDSP_Length(lhs.row), vDSP_Length(rhs.col), vDSP_Length(lhs.col))
        return Matrix(mul, row: lhs.row, col: rhs.col)
    }

    public func get() -> [[Float]] {
        var ret = [[Float]](repeating: [Float](repeating: 0.0, count: self.col), count: self.row)
        for i in 0..<self.row {
            for j in 0..<self.col {
                ret[i][j] = self.flat[i*self.col+j]
            }
        }
        return ret
    }
    public func transpose() -> Matrix {
        var trans = [Float](repeating: 0.0, count: self.row * self.col)
        vDSP_mtrans(self.flat, 1, &trans, 1, vDSP_Length(self.col), vDSP_Length(self.row))
        return Matrix(trans, row: self.col, col: self.row)
    }
    func compare(_ other:float4x4) -> Bool{
        return compare(other, 4,4)
    }
    func compare(_ other:float3x3) -> Bool{
        return compare(other, 3,3)
    }
    func compare(_ other:float2x2) -> Bool{
        return compare(other, 2,2)
    }
    private func compare(_ other: any SmallMatrix,_ row:Int,_ col:Int) -> Bool {
        guard self.row == row, self.col == col else {
            return false
        }
        for i in 0..<row {
            for j in 0..<col {
                let index = i * row + j
                if abs(self.flat[index] - other[i,j]) > 0.001 {
                    return false
                }
            }
        }
        return true
    }
}
protocol SmallMatrix {
    associatedtype Scalar: SIMDScalar
    subscript(row: Int, column: Int) -> Float { get set }
}
extension float4x4: SmallMatrix {
    typealias Scalar = Float
    subscript(_ row: Int,_  column: Int) -> Scalar {
        get { self[row][column] }
        set { self[row][column] = newValue }
    }
}
extension float3x3: SmallMatrix {
    typealias Scalar = Float
    subscript(_ row: Int,_  column: Int) -> Scalar {
        get { self[row][column] }
        set { self[row][column] = newValue }
    }
}
extension float2x2: SmallMatrix {
    typealias Scalar = Float
    subscript(_ row: Int,_  column: Int) -> Scalar {
        get { self[row][column] }
        set { self[row][column] = newValue }
    }
}
