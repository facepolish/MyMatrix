//
//  Complex.swift
//  MyMatrix
//
//  Created by 酒井裕司 on 2025/03/05.
//

import Foundation
protocol ComplexNumberType {
    associatedtype RealType: FloatingPoint //FloatingPoint
    var real: RealType { get set }
    var imaginary: RealType { get set }

    init(real: RealType, imaginary: RealType)
    func magnitude() -> RealType
    func argument() -> RealType
    static func +(lhs: Self, rhs: Self) -> Self
    static func -(lhs: Self, rhs: Self) -> Self
    static func *(lhs: Self, rhs: Self) -> Self
    static func /(lhs: Self, rhs: Self) -> Self
}

struct Complex<T: BinaryFloatingPoint>: ComplexNumberType {
    typealias RealType = T
    var real: T
    var imaginary: T

    init(real: T, imaginary: T) {
        self.real = real
        self.imaginary = imaginary
    }

    func magnitude() -> T {
        return (real * real + imaginary * imaginary).squareRoot()
    }

    func argument() -> T {
        return T(atan2(Double(imaginary), Double(real)))
    }

    static func +(lhs: Complex, rhs: Complex) -> Complex {
        return Complex(real: lhs.real + rhs.real, imaginary: lhs.imaginary + rhs.imaginary)
    }

    static func -(lhs: Complex, rhs: Complex) -> Complex {
        return Complex(real: lhs.real - rhs.real, imaginary: lhs.imaginary - rhs.imaginary)
    }

    static func *(lhs: Complex, rhs: Complex) -> Complex {
        return Complex(real: lhs.real * rhs.real - lhs.imaginary * rhs.imaginary, imaginary: lhs.real * rhs.imaginary + lhs.imaginary * rhs.real)
    }
    static func /(lhs: Complex, rhs: Complex) -> Complex {
        let denominator = rhs.real * rhs.real + rhs.imaginary * rhs.imaginary
        guard denominator != 0 else { fatalError("Division by zero") }
        return Complex(real: (lhs.real * rhs.real + lhs.imaginary * rhs.imaginary) / denominator, imaginary: (lhs.imaginary * rhs.real - lhs.real * rhs.imaginary) / denominator)
    }
}

