# Ran to time = 50000
dt = [6; 3; 2; 1.5];
%N = [1/4400; 2/4400; 3/4400; 4/4400; 8/4000];
N = [1/6423.2; 1/3103.5; 1/2046.3; 1/1526.4; 1/757.0248];

logN = log(N);


errH = [273.391120;54.039883; 22.990531; 12.900193; 3.158165];
errQx = [8171.200649;934.272114; 274.845045; 133.834016; 37.511973];
errQy = [8201.632332;979.581585; 281.932175; 132.881907; 34.963521];

logerrH = log(errH);
logerrQx = log(errQx);
logerrQy = log(errQy);
figure(1)
plot(logN, logerrH, 'o-');
title('Error in H');
figure(2);
plot(logN, logerrQx, 'o-');
title('Error in Qx');
figure(3)
plot(logN, logerrQy, 'o-')
title('Error in Qy')

slopeH = zeros(4,1);
slopeQx = zeros(4,1);
slopeQy = zeros(4,1);
for i = 1:4
    denom = logN(i+1) - logN(i);
    slopeH(i) = (logerrH(i+1) - logerrH(i))/denom;
    slopeQx(i) = (logerrQx(i+1) - logerrQx(i))/denom;
    slopeQy(i) = (logerrQy(i+1) - logerrQy(i))/denom;
end



