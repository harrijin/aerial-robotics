rng('shuffle');
errors = [];
for x = -1:0.1:1
    for y = -1:0.1:1
        for z = -1:0.1:1
            % Set RBI
            RBI_actual = euler2dcm([x*pi y*pi z*pi]);
            % Generate synthetic data
            vIMat = zeros(10,3);
            vBMat = zeros(10,3);
            for k = 1:10
                vIMat(k,:) = rand(1,3);
                vBMat(k,:) = (RBI_actual*vIMat(k,:)')';
            end
            % Find RBI using wahbaSolver and synthetic data
            RBI_wahba = wahbaSolver(ones(10,1), vIMat, vBMat);
            errors = [errors;norm(RBI_wahba-RBI_actual)];
        end
    end
end
max(errors)