package org.jenetics.example;

/**
 * Created by XjacX on 10/28/2015.
 */
public class ParameterTest {

    public static void main(String[] args) {

        System.out.println("parameter test");

        TSP2 tsp = new TSP2();
        RealFunction5Ackley ackley = new RealFunction5Ackley();
        RealFunction5Rastrigin rastrigin = new RealFunction5Rastrigin();

        //tsp.parameterTest(10);
        //ackley.parameterTest(10);
        rastrigin.parameterTest(10);



    }
}
