//
//  ComplexUnit.swift
//  MyMatrix
//
//  Created by 酒井裕司 on 2025/03/06.
//

import Testing
import Foundation
@testable import MyMatrix

struct ComplexTests {
    let epsilonF:Float = 0.0001 // 誤差を適切に設定して下さい
    let epsilonD:Double = 0.0001
    let epsilonF80:Float80 = 0.0001

    @Test func testArgumentForFloat() async {
        let complexNumber = Complex<Float>(real: 1.0, imaginary: 1.0)
        let result = complexNumber.argument()
        #expect( abs(result - Float.pi / 4) < epsilonF) // 45度
    }
    
    @Test func testArgumentForDouble() async{
        let complexNumber = Complex<Double>(real: 1.0, imaginary: 1.0)
        let result = complexNumber.argument()
        #expect(abs(result - Double.pi / 4) < epsilonD) // 45度
    }
    
    @Test func testArgumentForFloat80() async{
        let complexNumber = Complex<Float80>(real: 1.0, imaginary: 1.0)
        let result = complexNumber.argument()
        #expect(abs(result - Float80.pi / 4) < epsilonF80) // 45度
    }
    @Test func testArgumentForCGFloat() async{
        let complexNumber = Complex<CGFloat>(real: 1.0, imaginary: 1.0)
        let result = complexNumber.argument()
        #expect(Double(abs(result - CGFloat.pi / 4)) < epsilonD) // 45度
    }

    @Test func testArgumentForNegativeReal() async{
        let complexNumber = Complex<Double>(real: -1.0, imaginary: 1.0)
        let result = complexNumber.argument()
        #expect(abs(result - (3.0 * Double.pi) / 4) < epsilonD) // 135度
    }

    @Test func testArgumentForNegativeImaginary() async{
        let complexNumber = Complex<Float>(real: 1.0, imaginary: -1.0)
        let result = complexNumber.argument()
        #expect(abs(result + Float.pi / 4) < epsilonF) // -45度
    }

    @Test func testArgumentForNegativeRealAndImaginary() async{
        let complexNumber = Complex<Double>(real: -1.0, imaginary: -1.0)
        let result = complexNumber.argument()
        #expect(abs(result  + 3.0 * Double.pi / 4) < epsilonD) // -135度
    }
    
    @Test func testArgumentForPurelyImaginary() async{
        let complexNumber = Complex<Float>(real: 0.0, imaginary: 1.0)
        let result = complexNumber.argument()
        #expect(abs(result - Float.pi / 2) < epsilonF) // 90度
    }

    @Test func testArgumentForZero() async{
        let complexNumber = Complex<Double>(real: 0.0, imaginary: 0.0)
        let result = complexNumber.argument()
        #expect(abs(result) < epsilonD) // 定義により、引数は0
    }
}
