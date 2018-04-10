function s2plt(sp)
%S2PLT Plot a reflectance spectrum for sample 2
figure
semilogy(sp.spect(:,1),sp.spect(:,2),'.');
title(['Sample 2, \lambda = ' num2str(sp.lambda,3) 'nm']);
xlabel('\theta, deg');
ylabel('Reflectance');
end

