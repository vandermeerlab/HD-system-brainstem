function tests = TEST_quadrature
tests = functiontests(localfunctions);
end

function TestConvertQEStatesToAngle(testCase)

% generate fake state data
nSamples = 100;
state_list = randi(4, nSamples, 1); 
state_tvec = linspace(0, 5, nSamples);
state_tsd = tsd(state_tvec, state_list');

angle_tsd = ConvertQEStatesToAngle([], state_tsd);

% verify that first output sample is 0
verifyTrue(testCase, angle_tsd.data(1) == 0);

end



function setupOnce(testCase)

end

function teardownOnce(testCase)

end
