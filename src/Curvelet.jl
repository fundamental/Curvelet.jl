### This implementation is roughly in the order defined in
### `Uniform Discrete Curvelet Transform' by Truong T. Nguyen and Hervé Chauris
### Equations will be referenced out of this paper.

### Common Lisp comment style will be used here for clarity
module Curvelet
    export curveletTransform, inverseCurveletTransform
    using DSP


    ### Eqn. 13 defines β(t) s.t. it smoothly maps the range -1...1 to 0...1
    ### This function is used to specify the transition regions of the windows
    function _beta(t::AbstractFloat)
        ## While the absolute value is not included in any of the equations, it seems
        ## to be required to properly calculate values with t~ -1.0
        if(abs(t) <= 1)
            sqrt(abs(sum([-5/32 21/32 -35/32 35/32 1/2].*(t.^[7 5 3 1 0]))))
        elseif(1 <= t)
            1.0
        else
            0.0
        end
    end

    ### ηa & ηb are used to define the respective radial and angular falloff rates
    ### Within the paper the value of 0.15 was used for both, though recovery
    ### guarentees hold for more values.
    const ETA_A = 0.15
    const ETA_B = 0.15

    #### Windowing functions
    ###
    _w1p(t) = _beta((pi-abs(t))/(pi*ETA_A))
    _w0p(t) = _w1p(2*t*(1+ETA_A))
    _w0(w::Tuple{Float64,Float64}) = _w0p(w[1]).*_w0p(w[2])

    _w1(w::Tuple{Float64, Float64}) = ((1-_w0(w)^2)^0.5)*_w1p(w[1])*_w1p(w[2])

    _v1p(t::Float64,N) = _beta((2/N-1-t)/(2*ETA_B/N))*_beta((t+1)/(2ETA_B/N))

    _vp(t,N,l) = _v1p(t-(2*(l-1))/N,N)

    ## Eqn. 29 defines T(θ) as a warping function that warps the 1D function ṽ(t,l)
    ## to the angular function v(θ,l)
    function T(theta)
        if(-pi/4 <= theta && theta <= pi/4)
            tan(theta)
        elseif(pi/4 < theta && theta < pi/2)
            2-tan(pi/2-theta)
        elseif(-pi/2 < theta && theta < -pi/4)
            -2-tan(pi/2-theta)
        elseif(theta < -pi/2)
            -2
        else#  theta > pi/2
            2
        end
    end

    ### Warped version of ṽ(t,l)
    function _vl(theta,N,l)
        if(l<N)
            _vp(float(T(theta)), N, l)
        else
            _vp(T(pi/2-theta),N,2*N+1-l)
        end
    end

    ### Intermediate in defining the windowing function for l!=0
    function _ulp(w::Tuple{AbstractFloat, AbstractFloat}, N, l)
        _vl(atan2(w[2],w[1]),N,l)*_w1(w)
    end


    function _ul(w::Tuple{AbstractFloat, AbstractFloat}, N::Int, l::Int)
        if(l==0)
            _w0(w)
        else
            result = 0
            ## These extra loops should not be needed, but something is missed
            ## if these are ommited
            for a=[-2pi 0.0 2pi], b=[-2pi 0.0 2pi]
                result += _ulp((w[1]+a,w[2]+b), N, l)
            end
            result
        end
    end

    grid(x,size) = (-pi+2.0*pi/size*(x-1))

    ### Defines the window (without scaling term) exactly as described in the paper.
    ### This approach is quite expensive though, resulting ~100x the work really
    ### needed to get the output as many points outside the ROS are dealt with
    function direct_window_definition(size::Int, N::Int, l::Int)
        fn = (x,y)->_ul((x,y),N,l)

        result = zeros(size, size)

        ## Calculate the window, allowing it to have its full image expressed by
        ## wraping around the edge
        for x=1:size, y=1:size,
            result[y,x] = fn(grid(x,size), grid(y,size))
        end
        return result
    end

    ###Generate scaled window, see fast_window_definition in the optimization
    ### section
    function window(size::Int, N::Int, l::Int)
        if(l==0)
            fast_window_definition(size, N, l)*2
        else
            fast_window_definition(size, N, l)*2*sqrt(2)
        end
    end

    ### Performs the downsample by a factor of two in both directions
    function downsample(X)
        DSP.fftshift(DSP.fft(DSP.ifft(DSP.ifftshift(X))[1:2:end,1:2:end]))
    end

    function downsample2(X)
        Xx = DSP.fftshift(X)
        y  = X[1:Int(end/2),1:Int(end/2)]
        y += X[1:Int(end/2),Int(end/2)+1:end]
        y += X[Int(end/2)+1:end,Int(end/2)+1:end]
        y += X[Int(end/2)+1:end,1:Int(end/2)]
        DSP.fftshift(y)/4
    end

    function super_downsample(X)
        Xx = DSP.fftshift(X)
        y  = X[1:end/2,1:end/2]
        y += X[1:end/2,end/2+1:end]
        y += X[end/2+1:end,end/2+1:end]
        y += X[end/2+1:end,1:end/2]
        yy = y[1:end/2,:] + y[end/2+1:end,:]
        yyy = yy[1:end/2,:] + yy[end/2+1:end,:]
        DSP.fftshift(yyy)
    end

    function super_upsample(X)
        Xx = DSP.fftshift(X)
        XX = [Xx Xx; Xx Xx]
        XXX = [XX; XX]
        XXXX = [XXX; XXX]
    end


    function upsample(X)
        N = size(X,1)
        xx = zeros(2*N,2*N)+0im
        xx[1:2:end,1:2:end] = DSP.ifft(DSP.ifftshift(X))
        DSP.fftshift(DSP.fft(xx))
    end

    function upsample2(X)
        Xx = DSP.fftshift(X)
        [Xx Xx; Xx Xx]
    end


    C1 = 7
    C2 = 2*C1+1

    function _gen_windows(Wsize,N)
        win = Array{Matrix{Float64}}(2*N+1)
        for i=0:N
            print('.')
            win[i+1] = window(Wsize,N,i)
        end
        for i=N+1:2N
            win[i+1] = flipdim(win[i-N+1],1)'
            win[i+1][:,[2:end..., 1]] = win[i+1][:, [1:end...]]
        end
        win
    end

    struct CurveletPlan
        dim::Int
        angles::Int
        windows::Vector{Matrix{Float64}}
        subplan::Union{CurveletPlan,Void}
    end

    ## no need to print out a ton of matricies...
    import Base.show
    show(io::IO, cp::Curvelet.CurveletPlan) = print(io, "CurveletPlan")

    function planCurveletTransform(dim::Int, N::Int)
        print("<planning>")
        @assert mod(dim,2)==0
        @assert dim > 0
        @assert N > 0

        if(dim <= 16)
            return nothing
        end

        CurveletPlan(dim,2*N+1,_gen_windows(dim,N),planCurveletTransform(Int(dim/2), N))
    end

    struct CurveletCoeff
        HP::Vector{Matrix{Complex{Float64}}}
        LP::Union{CurveletCoeff,Matrix{Complex{Float64}}}
    end

    curveletTransform(x) = curveletTransform(x, planCurveletTransform(size(x,1), C1))

    function curveletTransform(x::Matrix{Float64}, plan::CurveletPlan)
        println("FFT")
        S = DSP.fftshift(DSP.fft(x))

        tf = Array(Matrix{Complex{Float64}}, plan.angles)
        print("Transforming")
        for i=1:plan.angles
            print(".")
            tf[i] = S.*plan.windows[i]
        end
        println()

        ds = Array(Matrix{Complex{Float64}}, plan.angles)
        print("Downsampling")
        for i=1:plan.angles
            ds[i] = downsample2(tf[i])
            print(".")
        end
        println()

        coeff = Array(Matrix{Complex{Float64}}, plan.angles)
        print("creating coeff")
        for i=1:plan.angles
            coeff[i] = DSP.ifft(DSP.fftshift(ds[i]))
            print(".")
        end
        println()


        CurveletCoeff(coeff[2:end], curveletTransform(ds[1], plan.subplan))
    end

    function curveletTransform(S::Matrix{Complex{Float64}}, plan::Void)
        DSP.ifft(DSP.fftshift(S))
    end

    function curveletTransform(S::Matrix{Complex{Float64}}, plan::CurveletPlan)
        tf = Array(Matrix{Complex{Float64}}, plan.angles)
        print("Transforming")
        for i=1:plan.angles
            print(".")
            tf[i] = S.*plan.windows[i]
        end
        println()

        ds = Array(Matrix{Complex{Float64}}, plan.angles)
        print("Downsampling")
        for i=1:plan.angles
            ds[i] = downsample2(tf[i])
            print(".")
        end
        println()

        coeff = Array(Matrix{Complex{Float64}}, plan.angles)
        print("creating coeff")
        for i=1:plan.angles
            coeff[i] = DSP.ifft(DSP.fftshift(ds[i]))
            print(".")
        end
        println()


        CurveletCoeff(coeff[2:end], curveletTransform(ds[1], plan.subplan))
    end

    inverseCurveletTransform(x) = inverseCurveletTransform(x, planCurveletTransform(size(x.HP[1],1)*2, C1))

    function inverseCurveletTransform(cc::CurveletCoeff, plan::CurveletPlan)
        Wsize = plan.dim

        us = Array{Matrix{Complex{Float64}}}(plan.angles)
        print("upsampling")
        for i=1:plan.angles-1
            us[i+1] = upsample2(DSP.fftshift(DSP.fft(cc.HP[i])))
            print(".")
        end
        println()
        us[1] = upsample2(inverseCurveletTransform(cc.LP, plan.subplan, true))
        println()

        prev = zeros(Wsize,Wsize)
        print("Inverting")
        for i=1:plan.angles
            print(".")
            prev += us[i].*plan.windows[i]
        end
        println()

        println("IFFT")
        println("imaginary : ",norm(imag(DSP.ifft(DSP.ifftshift(prev)))))
        real(DSP.ifft(DSP.ifftshift(prev)))
    end
    
    function inverseCurveletTransform(cc::CurveletCoeff, plan::CurveletPlan, asdf::Bool)
        Wsize = plan.dim

        us = Array{Matrix{Complex{Float64}}}(plan.angles)
        print("upsampling")
        for i=1:plan.angles-1
            us[i+1] = upsample2(DSP.fftshift(DSP.fft(cc.HP[i])))
            print(".")
        end
        println()
        #println("types are = ", typeof(cc.LP), " ", typeof(plan.subplan))
        us[1] = upsample2(inverseCurveletTransform(cc.LP, plan.subplan, true))
        println()

        prev = zeros(Wsize,Wsize)
        print("Inverting")
        for i=1:plan.angles
            print(".")
            prev += us[i].*plan.windows[i]
        end
        println()

        prev
    end

    function inverseCurveletTransform(cc::Matrix{Complex{Float64}}, plan::Void, asdf::Bool)
        DSP.fftshift(DSP.fft(cc))
    end

    #### Optimization methods
    ### While the overall transform as laid out above is quick, the generation of
    ### necessary windowing functions is not very fast if it is done with only first
    ### principle definitions.
    ### The overall windowing functions have a Region Of Support (ROS) which is
    ### mostly convex and slowly changing.
    ### This feature can be exploited to only calculate the points that are non-zero
    ### within the [-π,π]^2 square
    ###
    ### Through symmetry we can limit this optimization to the windows on the
    ### eastern side.
    ### From N&l we can estimate the approximate angle and use that to find a
    ### point within the ROS.
    ### From there the upper and lower bounds can be adjusted as the window is
    ### calculated column-wise.

    ### Fit the window bounds starting from the previous bounds
    function refit(f, low::Int, high::Int)
        while(f(low-1) != 0)
            low -= 1
        end
        while(f(low) == 0 && low != high)
            low += 1
        end
        while(f(high+1) != 0)
            high += 1
        end
        while(f(high) == 0 && high != low)
            high -= 1
        end

        return (low, high)
    end

    ## The oddness of the below code is either the result of some strange bug
    ## that I introduced or some WEIRD compiler bug with reguards to
    ## lambda precision...
    function run_cols(f, x, b, A, dir, N, l)
        #V = x::Float64->-pi+2.0*pi/size(A,1)*(x-1)
        M = size(A,1)
        (low, high) = refit(y->f(grid(x,M),grid(y,M)), b[1], b[2])

        for y=low:high
            x_::Int = mod(x-1, size(A,1))+1
            y_::Int = mod(y-1, size(A,2))+1
            A[y_, x_] = _ul((grid(x_,M),grid(y_,M)), N, l)
        end

        if(low != high)
            run_cols(f, x+dir, (low, high), A, dir, N, l)
        end
    end


    function fast_window_definition(size::Int, N::Int, l::Int)
        if(l>N)
            result = flipdim(fast_window_definition(size,N,N-l),1)'
            result[:,[2:end, 1]] = result[:, [1:end]]
            return result
        end
        fn = (x,y)->_ul((x,y),N,l)

        result = zeros(size, size)

        V = x->-pi+2*pi/size*(x-1)
        #f = (x,y)->fn(idxToVal(x), idxToVal(y))

        #-pi/2..3*pi/2 FIXME see eqn 28. T(theta)
        interp = x->(-1/2+1*x/N)
        slope  = (atan(pi*interp(l-1))+atan(pi*interp(l)))/2

        x::Int = l != 0 ? floor(size*4/5) : size/2 #TODO enhance this
        y::Int = round((size/(2*pi)*(V(x)*slope+pi))) + 1

        if(fn(V(x),V(y)) == 0)
            println("Warning: bad y value(",y,"), correcting...")
            ## Bad guess, which should not happen if the above math is
            ## right
            ## oh well, do a quick search to find a valid point
            for ny = 1:size
                if(fn(V(x),V(ny)) != 0)
                    y = ny
                    break
                end
            end
        end


        ## Real failure happens here
        @assert fn(V(x),V(y)) != 0

        bounds = refit(y->fn(V(x),V(y)), y, y)
        run_cols(fn, x, (y,y), result, -1, N, l)
        run_cols(fn, x+1, (y,y), result, +1, N, l)

        return result
    end
end
