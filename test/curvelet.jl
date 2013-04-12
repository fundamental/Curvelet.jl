using Images
peppers = imread("peppers.gif")
sample_data = reshape(peppers[1,:,:],512,512)'*1.0
##
###sample_data = rand(1024,1024)
win = windows(size(sample_data,1), C1)
###sample_data -= mean(sample_data)
####sample_data = data
@time tmp = transform(sample_data, win)
###ttmp = transform(ifft(ifftshift(tmp[1])))
####ttmp[1] = fftshift(fft(restore(transform(real(ifft(ifftshift(ttmp[1]))))))).*teh_window(32,3,0)
###println("mean is = ", mean(ttmp[1]))
###tmp_one_old = tmp[1]
####tmp[1] = fftshift(fft(restore(ttmp)))
@time sample_hat  = restore(tmp, win)
println("Done")
println("Error was: ", norm(sample_data-sample_hat))
