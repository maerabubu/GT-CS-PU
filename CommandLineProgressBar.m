classdef CommandLineProgressBar < handle
    %COMMANDLINEPROGRESSBAR
    % Progressbar for the command line, which also supports parfor loops.
    
    properties(Access=public)
        message = 'Progress'
        barLength = 42;
    end
    
    properties(Access=private)
        nIncrements_total
        nIncrements_current = 0
        wasAlreadyPrinted = false
        lastUpdate = -1
        dataQueue
    end
    
    properties(Dependent=true)
        progress
        progress_percent
    end
    
    methods
        function obj = CommandLineProgressBar(nIncrements_total)
            %COMMANDLINEPROGRESSBAR Construct an instance of this class
            
            % Validate input.
            validateattributes( ...
                nIncrements_total, ...
                {'numeric'}, ...
                {'real','finite','nonnan','nonsparse','nonempty','positive','integer','scalar'});
            
            % Assign property.
            obj.nIncrements_total = nIncrements_total;
            
            % Setup parallel pool.
            setupparallelpool
            
            % Create a dataqueue.
            obj.dataQueue = parallel.pool.DataQueue;
            
            % Set an after each event for the datacue.
            afterEach(obj.dataQueue, @(varargin) obj.update);
        end
        
        function increment(obj)
            % INCREMENT Trigger an update of the progressbar.
            send(obj.dataQueue,1);
        end
        
        function obj = update(obj)
            %UPDATE Update the progressbar.
            obj.nIncrements_current = obj.nIncrements_current+1;
            
            if obj.progress_percent ~= obj.lastUpdate
                obj.print;
                obj.lastUpdate = obj.progress_percent;
            end
        end
        
        function print(obj)
            %PRINT Print the current progress.
            if obj.barLength > 0
                nChars_complete = round((obj.barLength-2)*obj.progress);
                nChars_incomplete = obj.barLength-2-nChars_complete;
                
                barString = [ ...
                    '|' ...
                    repmat('=',1,nChars_complete) ...
                    repmat(' ',1,nChars_incomplete) ...
                    '| '
                    ];
            else
                barString = '';
            end
            
            % Print the message on the first increment. On all other
            % increments delte the old prograssbar chars.
            if obj.wasAlreadyPrinted
                fprintf(repmat('\b',1,numel(barString)+4));
            else
                fprintf('%s',obj.message)
                obj.wasAlreadyPrinted = true;
            end
            
            % Print the progressbar.
            fprintf('%s%3d%%',barString,obj.progress_percent);
            
            % On the last increment print a linebreak.
            if obj.progress == 1
                fprintf('\n');
            end
        end
        
        %% Getters
        function value = get.progress(obj)
            value = obj.nIncrements_current/obj.nIncrements_total;
        end
        
        function value = get.progress_percent(obj)
            value = floor(obj.progress*100);
        end
        
        %% Setters
        function set.message(obj,value)
            % Validate input.
            validateattributes( ...
                value, ...
                {'char'}, ...
                {'row','vector'});
            
            % Assign property.
            obj.message = value;
        end
        
        function set.barLength(obj,value)
            % Validate input.
            validateattributes( ...
                value, ...
                {'numeric'}, ...
                {'real','finite','nonnan','nonsparse','nonempty','integer','scalar','>',2});
            
            % Assign property.
            obj.barLength = value;
        end
    end
end

