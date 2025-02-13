//
//  Matrix.swift
//  Matrix
//
//  Created by 酒井裕司 on 2025/02/03.
//

import simd
import Foundation
import Accelerate

enum Errors: Error {
    case MatrixSizeMismatch
    case VectorSizeMismatch
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
        for i in 0..<self.col {
            let vec = self.vector(i, orientation: true)
            var difVec = Vector(0.0,size:vec.size)
            difVec = try vecs.reduce(difVec){try $0 - (try vec*$1) * $1}
            let retVec = try vec + difVec
            vecs.append(retVec.normalize())
        }
        return Matrix(vecs)
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
    
    public static func - (lhs: Matrix, rhs: Matrix) throws -> Matrix {
        guard lhs.row == rhs.row && lhs.col == rhs.col else{
            throw Errors.MatrixSizeMismatch
        }
        let minus = zip(lhs.flat, rhs.flat).map { $0.0 - $0.1 }
        return Matrix(minus, row: rhs.row, col: rhs.col)
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
