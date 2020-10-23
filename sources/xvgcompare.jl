# author : charlie
# date : 20201022

using UnicodePlots


function main()
    filenames = ARGS
    # println("filename -> $(filename)")

    datas = []
    legend_list_show = []
    xlabel_show = ""
    ylabel_show = ""
    column_num_show = 0
    
    for filename in filenames
        content = readlines(filename)
        content = [ line for line in content if length(line) != 0]
        column_num = length(split(content[end], ' ', keepempty=false))
    
        data = [[] for i in 1:column_num]
        legend_list = []
        xlabel = "default"
        ylabel = "default"
    
        for line in content
            if (line[1] in "@") == true
                if occursin("@    title", line)
                    fig_title = strip(strip(line[11:end], ' '), '"')
                    # println("fig_title -> $(fig_title)")
                elseif occursin("@    xaxis", line)
                    xlabel = strip(strip(line[18:end], ' '), '"')
                    # println("xlabel -> $(xlabel)")
                elseif occursin("@    yaxis", line)
                    ylabel = strip(strip(line[18:end], ' '), '"')
                    # println("ylabel -> $(ylabel)")
                elseif occursin("@ s", line) && occursin("legend", line)
                    push!(legend_list, strip(strip(line[12:end], ' '), '"'))
                    # println("legends -> $(legend_list)")
                end
            
            elseif (line[1] in "#@") == false
                line_data = split( line, ' ', keepempty = false)
                # for energy_compute.xvg
                if "time" in line_data      # && "ETOTAL" in line_data && "COULOMB" in line_data
                    legend_list = line_data[2:end]
                    xlabel = "--time(ps)--"
                    ylabel = "--(kJ/mol)--"
                    continue
                end
                # dump data into data variable
                for i in 1:column_num
                    push!(data[i], parse(Float64, line_data[i]))
                end
            end
        end

        push!(datas, data)
        if length(legend_list_show) == 0
            legend_list_show = legend_list
        end
        if legend_list_show != legend_list
            println(">>>>>>>>>>> legends of different file is not the same <<<<<<<<<")
        end
        if length(xlabel_show) == 0
            xlabel_show = xlabel
        end
        if xlabel_show != xlabel
            println(">>>>>>>>>>> xlabels of different file is not the same <<<<<<<<<")
        end
        if length(ylabel_show) == 0
            ylabel_show = ylabel
        end
        if ylabel_show != ylabel
            println(">>>>>>>>>>> ylabels of different file is not the same <<<<<<<<<")
        end
        if column_num_show == 0
            column_num_show = column_num
        end
        if column_num_show != column_num
            println(">>>>>>>>>>>>> WRONG, num of columns of different files is not the same <<<<<<<<<<<")
            println(">>>>>>> Exit <<<<<")
            exit(0)
        end
    end

    # diaw fig
    # in this comparation, legends should be the titles, filename should be the legend.
    if length(legend_list_show) != column_num_show-1
        legend_list_show = ["legend" for i in 2:column_num_show]
    end
    fig_width = 100
    for filename in filenames
        if length(ylabel_show*filename) > 10
            fig_width = 80
        end
    end
    for i in 2:column_num_show
        # guess the range of y 
        maxs = []
        mins = []
        for j in 1:length(filenames)
            push!(maxs, maximum(datas[j][i]))
            push!(mins, minimum(datas[j][i]))
        end
        ylim = [minimum(mins), maximum(maxs)]

        plt = lineplot([x for x in datas[1][1]], [y for y in datas[1][i]],
                       title = legend_list_show[i-1], name = filenames[1],
                       xlabel = xlabel_show, ylabel = ylabel_show, ylim = ylim,
                       width = fig_width, height = 20)
        for m in 2:length(filenames)
            lineplot!(plt, [x for x in datas[m][1]], [y for y in datas[m][i]], name = filenames[m])
        end
        println(plt)
    end
end

        





main()
