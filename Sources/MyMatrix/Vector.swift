//
//  Vector.swift
//  Matrix
//
//  Created by 酒井裕司 on 2025/02/03.
//

import Foundation
import Accelerate

public class Vector{
    var flat:[Float] = []
    var size:Int = 0
    init(_ a:[Float]){
        self.flat = a
        self.size = a.count
    }
    init(_ value:Float,size:Int){
        let valueFlat = [Float](repeating: value, count:size)
        self.flat = valueFlat
        self.size = size
    }
    func norm() -> Float {
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
    static func * (scalar:Float, vec:Vector) -> Vector {
        let ret = vec.flat.map {scalar * $0}
        return Vector(ret)
    }
    static func * (lhs:Vector,rhs:Vector)throws -> Float {
        guard lhs.size == rhs.size else{
            throw Errors.VectorSizeMismatch
        }
        var retValue:Float = 0
        vDSP_dotpr (lhs.flat,1,rhs.flat,1,&retValue,vDSP_Length(lhs.size))
        return retValue
    }
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

