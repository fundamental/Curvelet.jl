Curvelet.jl - The 2D Curvlet Transformation
==========================================
Mark McCurry
0.0.0, April 12th 2013


The curvelet transform is a fairly recent image processing technique that is
able to easily approximate curves present in images.
This package is an implementation of the ``Uniform Discrete Curvelet Transform''
as described in ``Uniform Discrete Curvelet Transform'' by Truong T. Nguyen and
Hervé Chauris.

Basic usage is as follows:

----------------------------------------
require("src/Curvelet.jl")

x = rand(128,128)
X = Curvelet.curveletTransform(x)
y = Curvelet.inverseCurveletTransform(X)
----------------------------------------

Restrictions
------------

Currently this transform works only for a simple class of inputs:
square images with dimensions that are powers of two in length and at least
16x16.
