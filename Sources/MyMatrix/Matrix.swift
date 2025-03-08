//
//  Matrix.swift
//  Matrix
//
//  Created by 酒井裕司 on 2025/02/03.
//

import simd
import Foundation
import Accelerate
// 定数近似値
let epsilon:Float = 0.0001
// エラーを定義する列挙体
enum Errors: Error {
    case MatrixSizeMismatch
    case VectorSizeMismatch
    case NoEigenValues
    case ParameterError
}
//typealias Scalar = Float
//typealias Scalar4x4 = float4x4

func eigenvaluesAndEigenvectors(matrix: [Double], order: Int) -> (eigenvalues: [(real: Double, imag: Double)], eigenvectors: [[Double]])? {
    // 行列が正方行列であることを前提に要素数との一致性を確認
    guard matrix.count == order * order else { return nil }
    // 固有値と固有ベクトルを格納するための配列
    var eigenvalueReal = [Double](repeating: 0.0, count: order)
    var eigenvalueImaginary = [Double](repeating: 0.0, count: order)
    var eigenvector = [Double](repeating: 0.0, count: order * order)
    var varMatrix = matrix
    var workspaceQuery = [Double](repeating: 0.0, count: 1)
    var status = Int32(0)
    // ワークスペースのクエリ
    var orderInt32 = Int32(order)
    var orderInt32_2 = Int32(order)

    var lda = Int32(order)
    var ldvr = Int32(order)
    var lwork = Int32(-1)
    
    // LAPACK関数を呼び出し、ワークスペースのサイズを求める
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
    // ステータスチェック
    guard status == 0 else { return nil }
    // 実際の計算
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
    // ステータスチェック
    guard status == 0 else { return nil }
    // 固有値と固有ベクトルの取得
    let eigenvalues = zip(eigenvalueReal, eigenvalueImaginary).map { (real: $0.0, imag: $0.1) }
    let eigenvectors = stride(from: 0, to: eigenvector.count, by: order).map { Array(eigenvector[$0..<$0+order]) }

    return (eigenvalues: eigenvalues, eigenvectors: eigenvectors)
}

public class Matrix<T:BinaryFloatingPoint> {
    // 行列のすべての要素を一つの配列にフラットに保存
    public var flat: [T] = []
    public var row: Int // 行数
    public var col: Int // 列数
    
    // 初期化子：2次元配列を受け取りMatrixを生成
    public init(_ a: [[T]]) throws {
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
    // 行列の指定位を取得
    public func get(_ row:Int,_ col:Int) -> T {
        let offset = self.col * row + col
        return self.flat[offset]
    }

    public init(_ a: [T], row: Int, col: Int) {
        self.flat = a
        self.row = row
        self.col = col
    }
    public init(_ a:Vector<T>,orientation:Bool){//if orientation is true -> column vector
        self.flat = a.flat
        if orientation {
            self.row = a.size
            self.col = 1
        }else{
            self.row = 1
            self.col = a.size
        }
    }

    public init(_ f:[Vector<T>]){// assume vecgors as col vectors
        let cols = f.count
        let rows = f[0].size
        let flat_ini:[T] = []
        let flat = f.reduce(flat_ini){$0+$1.flat}
        let retMat = Matrix(flat,row:cols,col:rows).transpose()
        self.flat = retMat.flat
        self.row = retMat.row
        self.col = retMat.col
    }
    public func orthnormal() throws -> Matrix {
        var vecs:[Vector<T>] = []
        let count = min(self.col,self.row)
        for i in 0..<count {
            let vec = self.vector(i, orientation: true)
            var difVec = Vector(T(0.0),size:vec.size)
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
    // 特定のインデックスの列ベクトルまたは行ベクトルを取得
    //extract row or column vector from matrix. if orientation is true -> column vec
    public func vector(_ index:Int,orientation:Bool) -> Vector<T> {
        var flat:[T] = []
        if orientation {// Column vector
            flat = stride(from:index,to:self.flat.count,by:self.col).map{i in self.flat[i]}
        }else{
            let offset = index*self.col
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
    // 非負行列可動かを判定
    public func isNonNegativeMatrix() -> Bool {
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
    //対角成分の抽出
    public func mainDiagonal() -> [T] {
        var ret:[T] = []
        let count = min(self.row,self.col)
        for i in 0..<count {
            ret.append(self.get(i,i))
        }
        return ret
    }
    public func googleVector() throws -> (eigen:T,eigenVec:Vector<T>)? {// If the matrix is no-negative google vector exist
        guard self.row == self.col else {
            print ("Matrix is not square")
            throw Errors.ParameterError
        }
        guard self.isNonNegativeMatrix() else {
            print ("Matrix is not Non-Negative")
            throw Errors.ParameterError
        }
        let eigens = try self.realEigen()
        guard eigens.eigenVectors.count != 0 else {
            print ("Matrix has no eigenvalue")
            return nil
        }
        var maxVal:T = 0
        var maxVec:Vector<T> = eigens.eigenVectors[0]
        for i in 0..<eigens.eigenValues.count {
            if eigens.eigenValues[i] > maxVal {
                maxVal = eigens.eigenValues[i]
                maxVec = eigens.eigenVectors[i]
            }
            
        }
        return (maxVal,maxVec)
    }
    //実数固有値の抽出
    public func realEigen() throws -> (eigenValues:[T], eigenVectors:[Vector<T>]){// calculate eigen values for real matrix
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
        var retVectors:[Vector<T>] = []
        var eigens:[T] = []
        
        for item in zip(result.eigenvalues,result.eigenvectors){
            if item.0.imag != 0 {// ignore imaginary eigen values
//                eigens.append(0)
            }else {
                eigens.append(T(item.0.real))
                let vec = Vector(item.1.map{T($0)})
                retVectors.append(vec)
            }
        }
        return (eigens,retVectors)
    }
    public static func - (lhs: Matrix, rhs: Matrix) throws -> Matrix {
        guard lhs.row == rhs.row && lhs.col == rhs.col else{
            throw Errors.MatrixSizeMismatch
        }
        let minus = zip(lhs.flat, rhs.flat).map { $0.0 - $0.1 }
        return Matrix(minus, row: rhs.row, col: rhs.col)
    }
    public static func * (lhs:Matrix, rhs:Vector<T>) throws -> Vector<T> { // Assume vector as col vector
        guard lhs.col == rhs.size else {
            throw Errors.MatrixSizeMismatch
        }
        let matt = lhs.get()
        var flat:[T] = []
        for i in 0..<lhs.row {
            var val:T = 0
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
        var mul = [T](repeating: 0.0, count: lhs.row * rhs.col)
        if T.self == Float.self {
            let lhsFloat = lhs.flat as! [Float]
            let rhsFloat = rhs.flat as! [Float]
            var resultFloat = [Float](repeating: 0.0, count: mul.count)
            vDSP_mmul(lhsFloat, 1, rhsFloat, 1, &resultFloat, 1, vDSP_Length(lhs.row), vDSP_Length(rhs.col), vDSP_Length(lhs.col))
            mul = resultFloat as! [T]
        } else if T.self == Double.self {
            let lhsFloat = lhs.flat as! [Double]
            let rhsFloat = rhs.flat as! [Double]
            var resultFloat = [Double](repeating: 0.0, count: mul.count)
            vDSP_mmulD(lhsFloat, 1, rhsFloat, 1, &resultFloat, 1, vDSP_Length(lhs.row), vDSP_Length(rhs.col), vDSP_Length(lhs.col))
            mul = resultFloat as! [T]
        }
        else {
                   fatalError("Unsupported floating-point type")
               }
        return Matrix(mul, row: lhs.row, col: rhs.col)
    }

    public func get() -> [[T]] {
        var ret = [[T]](repeating: [T](repeating: 0.0, count: self.col), count: self.row)
        for i in 0..<self.row {
            for j in 0..<self.col {
                ret[i][j] = self.flat[i*self.col+j]
            }
        }
        return ret
    }
    //転置演算
    public func transpose() -> Matrix {
        var trans = [T](repeating: 0.0, count: self.row * self.col)
        if T.self == Float.self {
            let flat = self.flat as! [Float]
            var transFloat = [Float](repeating: 0.0, count: self.row * self.col)
            vDSP_mtrans(flat, 1, &transFloat, 1, vDSP_Length(self.col), vDSP_Length(self.row))
            trans = transFloat as! [T]
        }else if T.self == Double.self {
            let flat = self.flat as! [Double]
            var transFloat = [Double](repeating: 0.0, count: self.row * self.col)
            vDSP_mtransD(flat, 1, &transFloat, 1, vDSP_Length(self.col), vDSP_Length(self.row))
            trans = transFloat as! [T]
        }
        return Matrix(trans, row: self.col, col: self.row)
    }
}
func multiplyMatrix(_ a:Matrix<Float>,_ b:Matrix<Float>) -> Matrix<Float>{
    let a_flat = a.flat.map{Double($0)}
    let b_flat = b.flat.map{Double($0)}
    let ret = DGEMM(a: a_flat, b: b_flat, m: a.row, n: b.col, k: a.col)
    let fl = ret.map{Float($0)}
    return Matrix(fl,row:a.row,col:b.col)
}
func multiplyMatrix(_ a:Matrix<Double>,_ b:Matrix<Double>) -> Matrix<Double>{
    let ret = DGEMM(a: a.flat, b: b.flat, m: a.row, n: b.col, k: a.col)
    return Matrix(ret,row:a.row,col:b.col)
}
func DGEMM(
    a: [Double],
    b: [Double],
    m: Int,
    n: Int,
    k: Int
) -> [Double] {
    let order = CblasRowMajor
    let trans = CblasNoTrans
    
    var c = [Double](repeating: 0.0, count: m * n)
    
    a.withUnsafeBufferPointer { aPtr in
        b.withUnsafeBufferPointer { bPtr in
            c.withUnsafeMutableBufferPointer { cPtr in
                cblas_dgemm(
                    order,
                    trans,
                    trans,
                    Int32(m),
                    Int32(n),
                    Int32(k),
                    1.0,
                    aPtr.baseAddress!,
                    Int32(k),
                    bPtr.baseAddress!,
                    Int32(n),
                    0.0,
                    cPtr.baseAddress!,
                    Int32(n)
                )
            }
        }
    }
    return c
}
// 行列の逆行列を計算する関数
func invert(matrix : [Double]) -> [Double] {
    var inMatrix = matrix
    // 行列のサイズを取得し、ローカル変数にコピー
    var N = __CLPK_integer(sqrt(Double(matrix.count)))
    var pivots = [__CLPK_integer](repeating: 0, count: Int(N))
    // 作業領域としてのバッファを用意
    var workspace = [Double](repeating: 0.0, count: Int(N))
    var error : __CLPK_integer = 0

    withUnsafeMutablePointer(to: &N) {
        // LU分解を実行
        dgetrf_($0, $0, &inMatrix, $0, &pivots, &error)
        // 逆行列を計算
        dgetri_($0, &inMatrix, $0, &pivots, &workspace, $0, &error)
    }
    return inMatrix
}
