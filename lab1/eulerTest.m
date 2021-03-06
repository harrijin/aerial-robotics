% test that euler angles can be converted to R and back
for i = -0.49:0.05:0.49
    for j = -1:0.1:1
        for k = -1:0.1:1
            % test that euler angles can be converted to R and back
            e = [pi*i;pi*j;pi*k];
            R = euler2dcm(e);
            e2 = dcm2euler(R);
            if norm(e2-e) > 0.0001
                e
                e2
                error("Did not convert properly");
            end
            % test that rotationMatrix(e) == R
            R2 = rotationMatrix([0;1;0], e(2))*rotationMatrix([1;0;0], e(1))*rotationMatrix([0;0;1], e(3));
            if norm(R2-R) > 0.0001
                error("Did not get same transformation");
            end
        end
    end
end