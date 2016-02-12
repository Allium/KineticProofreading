# KineticProofreading

Code for kinetic proofreading and sorting project with @DanSeeto.

## Purpose

The thermodynamics of "sorting" is an interesting topic from the perspective of fundamental physics, as well as its e.g. biological applications. In this project we simulate a process which attempts to sort particles of two species into two different boxes.

## The Setup

Two boxes connected by two channels. Each channel has a sorting device, which is tuned to permit passage of one type of particle in one direction only.

## The Sorting Devices

The sorting device we use is a Hopfield-style kinetic proofreader, i.e. we simulate the Hopfield reaction network. "Product" particles are transferred to the other box.

The sorting process consumes energy, and this is one thing we are interested in measuring.

## The Particles

The two species of particles differ in their unbinding rates by a factor &Delta; = Exp[-&Delta; F/T].

## More coming...