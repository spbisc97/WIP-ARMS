function plot_xu(x, u, time_array, names, xd,defects,order, type,plot_duration,plot_figure)

            %plot_xu(x, u, time_array, names,xd,order)
            if type == ""
                return
            end

            if nargin < 3
                error("missing entries")
            end

            [h, ~] = size(x);

            if nargin < 4 || isempty(names) || length(names)<h
                init=length(names)+1;
                names =[names,init:h]; %repmat("", 1, h);
            end


            desired = true;

            if nargin < 5
                desired = false;
            end

            if nargin < 7||isempty(order)

                order = (1:h)';
            end

            if nargin < 8
                type = 'GNMS_iLQR';
            end
            if nargin < 9
                plot_duration = inf;
            end
            if nargin < 10 || isempty(plot_figure)% || class(plot_figure)~="matlab.ui.Figure"
                figure_name="GNMS_iLQR";
            else
                set(0,'currentfigure',plot_figure)
                figure_name=plot_figure.Name;
            end

            
            [h, l] = size(order);

            if isempty(defects)
                tiles_n=h+1;
                def=false;
            else
                tiles_n=h+2;
                def=true;
                defplot  = double(any(circshift(defects,1,2)~=0,1));
                defplot(defplot==0)=nan;
            end





            tiledlayout(tiles_n, 1,'Padding','Compact','Tilespacing','Compact');

            for i = 1:h
                nexttile(i)
                j = 1;
                lgd = [];

                while j <= l && ~isnan(order(i, j))

                    idx = order(i, j);
                    plot(time_array, x(idx, :))
                    hold on
                    lgd = [lgd, names(idx)+""]; %#ok

                    if def
                        plot(time_array, x(idx, :).*defplot,"--o")
                        lgd = [lgd, ""]; %#ok
                    end
                    if desired && length(xd(:, 1)) >= idx
                        plot(time_array, xd(idx, :), LineStyle = "--")
                        hold on
                        lgd = [lgd, names(idx) + "_d"]; %#ok
                    end

                    title(type)

                    j =j+1;




                end
                legend(lgd,'Location','northeastoutside'); %interpreter =latex
                ylim('padded')
                xlim([-inf time_array(end)]);
                grid minor


            end
            xl = xlim;
            lgd=[];
            if def
                nexttile(h+1)
                plot(time_array(1:end), defects)
                for i=names
                    lgd=[lgd,"def-"+i];
                end
                legend(lgd,'Location','northeastoutside')
                title(type+"  defects")
                ylim('padded')
                xlim(xl);
                grid minor
            end




            nexttile(tiles_n)
            plot(time_array(1:end - 1), u)
            %legend("controls",'Location','northeastoutside')
            title(type+"   controls")
            ylim('padded')
            xlim(xl);
            grid minor

            if isinf(plot_duration)
                pause
            else
                pause(plot_duration)
            end
            set(gcf, 'Name', figure_name, 'NumberTitle', 'off')
        end %function