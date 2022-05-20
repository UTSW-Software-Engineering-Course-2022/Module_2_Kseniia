%Example
classdef Test1 < matlab.unittest.TestCase
    properties
        u
        m
        r
        l
        Fi
    end
    
    methods(TestMethodSetup)
        function CreateVars(testCase)
            u=ComUnit('erg',ComUnit.nm_to_cm(1000),300,ComUnit.kBT_to_erg(10,300)); 
            m=ModMembrane(2,'unit',u);
            r = mean(sqrt(sum(m.var.coord(:,1).^2+m.var.coord(:,2).^2+m.var.coord(:,3).^2,2)));
            [Fi] = Finternal(m, 'plot_or_not', true)
            
        
            testCase.u = u;
            testCase.m = m;
            testCase.Area_r = 4*pi*r^2
            testCase.l = edge_length(m.var.coord, m.var.edge_all);
            testCase.Fi = Fi
            testCase.X_idx, testCase.f_of_l = X_idx_and_f_of_l(Fi, l);
            
        end
    end
    
    methods (Test)
        function testEdgeLegth(testCase)
            actSolution = sum(edge_length(testCase.m.var.coord, testCase.m.var.edge_all));
            expSolution = 0
            testCase.verifyGreaterThanOrEqual(actSolution, expSolution)
        end
        
        function testAreaCompare(testCase)
            Area1 = testCase.m.Area()
            testCase.assertEqual(round(sum(Area1)), round(testCase.Area_r));
        end
        
        function testDeltaTFinal(testCase)
            Ftotal = comp_Ftotal(testCase.m, testCase.m.var.coord, testCase.f_of_l, testCase.l);
            delta_t_final = comp_delta_final(testCase.X_idx, testCase.m, testCase.l, testCase.Fi, testCase.m.var.coord, Ftotal);
            deltaThreshold = 0.1
            testCase.verifyLessThan(delta_t_final, deltaThreshold);
        end
    end
    
end

