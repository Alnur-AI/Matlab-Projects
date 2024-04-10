function solve = solve_problem(inputTaskTableData, conditionsTableData, segBegin, segEnd, stepsCount, accuracyExternal, accuracyInternal, solvingMethodForInternalTask, internalTaskStepSize, solvingMethodForExternalTask, timeOfT, initialVectorTableData, isStep)
        ftx = inputTaskTableData;
        ftx = ftx(1:end, 2);
        n = size(ftx, 1);
        
        
        % Находим производную
        symArray = sym('x%d', [1 n]);
        syms(symArray);
        tmpArr{1, n} = [];
        for i = 1:n
            tmpArr = [tmpArr, str2sym(ftx(i))];
        end
        dftx = jacobian(tmpArr, symArray);
        
        % Генерируем матрицу X
        tmpSymArray = zeros(0);
        XMatrix = zeros(0);
        for i = 1:n
            for j = 1:n
                tmpStr = ['x' num2str(n * i + j)];
                tmpSymArray = [tmpSymArray; sym(tmpStr)];
            end
            XMatrix = [XMatrix tmpSymArray];
            tmpSymArray = zeros(0);
        end
        
        % Генерируем начальные значения для матрицы Х
        initConditionForXMatrix = reshape(eye(n, n), [n*n, 1])';
        % Объединяем вектор начальных значений и начальные значения для матрицы Х
        initConditionsForInternalTask = [initialVectorTableData initConditionForXMatrix];
        % Формируем внутреннюю задачу
        DX = dftx*XMatrix;
        strDX = zeros(0);
        tmpStrDX = zeros(0);
        
        for i = 1:n
            for j = 1:n
                a = char(DX(j, i));
                A = string(a);
                tmpStrDX = [tmpStrDX; A];
            end
            strDX = [strDX tmpStrDX];
            tmpStrDX = zeros(0);
        end
        strDX = reshape(strDX, [n*n, 1]);
        toFile = [ftx; strDX];
        % Записываем задачу в файл
        save_sys_to_file(toFile);
        
        
        charFP = conditionsTableData;
        % Получаем Ф(р)
        m = size(charFP, 1);
        FP = zeros(0);
        tmpFP = zeros(0);

        for i = 1:m
            strTmp = charFP(i, 1);
            for j = 1:m
                toReplace = ['x' num2str(m + 1 - j) '(' num2str(segBegin) ')'];
                replacement = ['xa' num2str(m + 1 - j)];
                strTmp = strrep(strTmp, toReplace, replacement); 

                toReplace = ['x' num2str(m + 1 - j) '(' num2str(segEnd) ')'];
                replacement = ['xb' num2str(m + 1 - j)];
                strTmp = strrep(strTmp, toReplace, replacement); 

                toReplace = ['x' num2str(m + 1 - j) '(a)'];
                replacement = ['xa' num2str(m + 1 - j)];
                strTmp = strrep(strTmp, toReplace, replacement); 

                toReplace = ['x' num2str(m + 1 - j) '(b)'];
                replacement = ['xb' num2str(m + 1 - j)];
                strTmp = strrep(strTmp, toReplace, replacement); 
            end
            tmpFP = [tmpFP; string(strTmp)];
        end

        symArray = sym('xa%d', [1 m]);
        syms(symArray);
        symArray = sym('xb%d', [1 m]);
        syms(symArray);

        [token, remain] = strtok(tmpFP, '=');
        token = char(token);
        remain = char(remain);
        remain = remain(:, 2:end);

        for i = 1:m
            %FP = [FP; sym(token(i, :)) - sym(remain(i, :))];
            FP = [FP; str2sym(token(i, :)) - str2sym(remain(i, :))];
        end

        
        p0 = initialVectorTableData;
        
        global DRXA XAP DRXB XBP FP0;
        symArray1 = sym('xa%d', [1 n]);
        symArray1 = [symArray1 sym('xb%d', [1 n])];
        syms(symArray1);
        symArray1 = num2cell(symArray1);
        
        symArray2 = sym('xa%d', [1 n]);
        syms(symArray2);
        
        symArray3 = sym('xb%d', [1 n]);
        syms(symArray3);
        
        stepSize = 1;
        
        if isStep == 0
           tmp1 = [segBegin, segEnd];
           tmp2 = [segEnd, segBegin];
           tmp3 = [timeOfT, segBegin];
           tmp4 = [timeOfT, segEnd];
        else
           tmp1 = [segBegin:internalTaskStepSize:segEnd];
           tmp2 = [segEnd:internalTaskStepSize:segBegin];
           tmp3 = [timeOfT:-internalTaskStepSize:segBegin];
           tmp4 = [timeOfT:internalTaskStepSize:segEnd];
        end
        
        tStart = tic;
        waiting = waitbar(0,'Вычисляем...');
        for i = 1:stepsCount
            for j = 0:stepSize:(1-stepSize)
                switch solvingMethodForInternalTask
                    % Выбираем метод решения для внутренней задачи
                    case 1
                        if timeOfT == segBegin
                            [T, X] = ode45(@systemTemp, tmp1, initConditionsForInternalTask, odeset('RelTol', accuracyInternal));
                        elseif timeOfT == segEnd
                            [T, X] = ode45(@systemTemp, tmp2, initConditionsForInternalTask, odeset('RelTol', accuracyInternal));
                            T = T(end:-1:1, :);
                            X = X(end:-1:1, :);
                        else 
                            [T1, X1] = ode45(@systemTemp, tmp3, initConditionsForInternalTask, odeset('RelTol', accuracyInternal));
                            [T2, X2] = ode45(@systemTemp, tmp4, initConditionsForInternalTask, odeset('RelTol', accuracyInternal));
                            T = [T1(end:-1:2); T2];
                            X = [X1(end:-1:2, :); X2];
                        end
                    case 2
                        if timeOfT == segBegin
                            [T, X] = ode23(@systemTemp, tmp1, initConditionsForInternalTask, odeset('RelTol', accuracyInternal));
                        elseif timeOfT == segEnd
                            [T, X] = ode23(@systemTemp, tmp2, initConditionsForInternalTask, odeset('RelTol', accuracyInternal));
                            T = T(end:-1:1, :);
                            X = X(end:-1:1, :);
                        else 
                            [T1, X1] = ode23(@systemTemp, tmp3, initConditionsForInternalTask, odeset('RelTol', accuracyInternal));
                            [T2, X2] = ode23(@systemTemp, tmp4, initConditionsForInternalTask, odeset('RelTol', accuracyInternal));
                            T = [T1(end:-1:2); T2];
                            X = [X1(end:-1:2, :); X2];
                        end
                    case 3
                        if timeOfT == segBegin
                            [T, X] = ode23t(@systemTemp, tmp1, initConditionsForInternalTask, odeset('RelTol', accuracyInternal));
                        elseif timeOfT == segEnd
                            [T, X] = ode23t(@systemTemp, tmp2, initConditionsForInternalTask, odeset('RelTol', accuracyInternal));
                            T = T(end:-1:1, :);
                            X = X(end:-1:1, :);
                        else 
                            [T1, X1] = ode23t(@systemTemp, tmp3, initConditionsForInternalTask, odeset('RelTol', accuracyInternal));
                            [T2, X2] = ode23t(@systemTemp, tmp4, initConditionsForInternalTask, odeset('RelTol', accuracyInternal));
                            T = [T1(end:-1:2); T2];
                            X = [X1(end:-1:2, :); X2];
                        end
                    case 4
                        if timeOfT == segBegin
                            [T, X] = ode113(@systemTemp, tmp1, initConditionsForInternalTask, odeset('RelTol', accuracyInternal));
                        elseif timeOfT == segEnd
                            [T, X] = ode113(@systemTemp, tmp2, initConditionsForInternalTask, odeset('RelTol', accuracyInternal));
                            T = T(end:-1:1, :);
                            X = X(end:-1:1, :);
                        else 
                            [T1, X1] = ode113(@systemTemp, tmp3, initConditionsForInternalTask, odeset('RelTol', accuracyInternal));
                            [T2, X2] = ode113(@systemTemp, tmp4, initConditionsForInternalTask, odeset('RelTol', accuracyInternal));
                            T = [T1(end:-1:2); T2];
                            X = [X1(end:-1:2, :); X2];
                        end
                    case 5
                        if timeOfT == segBegin
                            [T, X] = ode23s(@systemTemp, tmp1, initConditionsForInternalTask, odeset('RelTol', accuracyInternal));
                        elseif timeOfT == segEnd
                            [T, X] = ode23s(@systemTemp, tmp2, initConditionsForInternalTask, odeset('RelTol', accuracyInternal));
                            T = T(end:-1:1, :);
                            X = X(end:-1:1, :);
                        else 
                            [T1, X1] = ode23s(@systemTemp, tmp3, initConditionsForInternalTask, odeset('RelTol', accuracyInternal));
                            [T2, X2] = ode23s(@systemTemp, tmp4, initConditionsForInternalTask, odeset('RelTol', accuracyInternal));
                            T = [T1(end:-1:2); T2];
                            X = [X1(end:-1:2, :); X2];
                        end
                    case 6
                        if timeOfT == segBegin
                            [T, X] = ode15s(@systemTemp, tmp1, initConditionsForInternalTask, odeset('RelTol', accuracyInternal));
                        elseif timeOfT == segEnd
                            [T, X] = ode15s(@systemTemp, tmp2, initConditionsForInternalTask, odeset('RelTol', accuracyInternal));
                            T = T(end:-1:1, :);
                            X = X(end:-1:1, :);
                        else 
                            [T1, X1] = ode15s(@systemTemp, tmp3, initConditionsForInternalTask, odeset('RelTol', accuracyInternal));
                            [T2, X2] = ode15s(@systemTemp, tmp4, initConditionsForInternalTask, odeset('RelTol', accuracyInternal));
                            T = [T1(end:-1:2); T2];
                            X = [X1(end:-1:2, :); X2];
                        end
                end

                n = size(FP, 1);
                
                
                % Находим Х(а, р)
                XAP = X(1, n+1:end);
                XAP = reshape(XAP, [n n]);
                
                % Находим Х(b, р)
                XBP = X(end, n+1:end);
                XBP = reshape(XBP, [n n]);

                XP0 = num2cell([X(1, 1:n) X(end, 1:n)]);

                FP0 = subs(FP(:, 1), symArray1, XP0);
                
                
                % Находим R'x
                DRXA = jacobian(FP, symArray2);
                % Находим R'y
                DRXB = jacobian(FP, symArray3);
                
                DRXA = subs(DRXA, symArray1, XP0);
                DRXB = subs(DRXB, symArray1, XP0);
                
                
                switch solvingMethodForExternalTask
                    % Выбираем метод решения для внешней задачи
                    case 1
                        [mu, p] = ode45(@findMatrixFPrime, [j (j+stepSize)], p0, odeset('RelTol', accuracyExternal));
                    case 2
                        [mu, p] = ode23(@findMatrixFPrime, [j (j+stepSize)], p0, odeset('RelTol', accuracyExternal));
                    case 3
                        [mu, p] = ode23t(@findMatrixFPrime, [j (j+stepSize)], p0, odeset('RelTol', accuracyExternal));
                    case 4
                        [mu, p] = ode113(@findMatrixFPrime, [j (j+stepSize)], p0, odeset('RelTol', accuracyExternal));
                    case 5
                        [mu, p] = ode23s(@findMatrixFPrime, [j (j+stepSize)], p0, odeset('RelTol', accuracyExternal));
                    case 6
                        [mu, p] = ode15s(@findMatrixFPrime, [j (j+stepSize)], p0, odeset('RelTol', accuracyExternal));
                end
                
                p(end, :);
                p0 = p(end, :);
                initConditionsForInternalTask = [p0 initConditionForXMatrix];
                waitbar((j / stepsCount) + ((i - 1) / stepsCount));
               
            end
        end
        close(waiting);
        tEnd = toc(tStart);
        % Получаем результат
        solve = [T(:, :), X(:, 1:n)];
end