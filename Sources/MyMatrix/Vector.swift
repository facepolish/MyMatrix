//
//  Vector.swift
//  Matrix
//
//  Created by 酒井裕司 on 2025/02/03.
//

import Foundation
import Accelerate

public struct Vector<T:BinaryFloatingPoint>{
    var flat:[T] = []
    var size:Int = 0
    init(_ a:[T]){
        self.flat = a
        self.size = a.count
    }
    init(_ value:T,size:Int){
        let valueFlat = [T](repeating: value, count:size)
        self.flat = valueFlat
        self.size = size
    }
    func norm() -> T {
        let ret = self.flat.reduce(0){$0+$1*$1}
        return sqrt(ret)
    }
    func normalize() -> Vector {
        let norm = self.norm()
        if norm == 0 {
            return self
        }
        let coef = 1/norm
        return coef * self
    }
    public func get() -> [T] {
        return self.flat
    }
    static func * (scalar:T, vec:Vector) -> Vector {
        let ret = vec.flat.map {scalar * $0}
        return Vector(ret)
    }
/*
    static func * (lhs:Vector,rhs:Vector)throws -> T {
        guard lhs.size == rhs.size else{
            throw Errors.VectorSizeMismatch
        }
        var retValue:T = 0
        vDSP_dotpr (lhs.flat,1,rhs.flat,1,&retValue,vDSP_Length(lhs.size))
        return retValue
    }
 */
    static func + (lhs:Vector,rhs:Vector)throws -> Vector {
        guard lhs.size == rhs.size else{
            throw Errors.VectorSizeMismatch
        }
        let retFlat = zip(lhs.flat,rhs.flat).map{$0.0 + $0.1}
        return Vector(retFlat)
    }
    static func - (lhs:Vector,rhs:Vector)throws -> Vector {
        guard lhs.size == rhs.size else{
            throw Errors.VectorSizeMismatch
        }
        let retFlat = zip(lhs.flat,rhs.flat).map{$0.0 - $0.1}
        return Vector(retFlat)
    }
}
extension Vector where T:BinaryFloatingPoint {
    static func * (lhs: Vector, rhs: Vector) throws -> T {
        guard lhs.size == rhs.size else {
            throw Errors.VectorSizeMismatch
        }
        // 型チェックと分岐
        if T.self == Float.self {
            let lhsFloat = lhs.flat as! [Float]
            let rhsFloat = rhs.flat as! [Float]
            var result: Float = 0
            vDSP_dotpr(
                lhsFloat, 1,
                rhsFloat, 1,
                &result,
                vDSP_Length(lhs.size)
            )
            return T(result)
            
        } else if T.self == Double.self {
            let lhsDouble = lhs.flat as! [Double]
            let rhsDouble = rhs.flat as! [Double]
            var result: Double = 0
            vDSP_dotprD(
                lhsDouble, 1,
                rhsDouble, 1,
                &result,
                vDSP_Length(lhs.size)
            )
            return T(result)
        } else {
            fatalError("Unsupported floating-point type")
        }
    }
}
