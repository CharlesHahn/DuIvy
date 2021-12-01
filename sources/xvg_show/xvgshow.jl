# author : charlie
# date : 20201022

using UnicodePlots


function main()
    filename = ARGS[1]
    # println("filename -> $(filename)")

    content = readlines(filename)
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
            if "time" in line_data && "ETOTAL" in line_data && "COULOMB" in line_data
                legend_list = line_data[2:end]
                xlabel = "time(ps)"
                ylabel = "(kJ/mol)"
                continue
            end
            # dump data into data variable
            for i in 1:column_num
                push!(data[i], parse(Float64, line_data[i]))
            end
        end
    end

    if length(legend_list) != column_num-1
        legend_list = ["legend" for i in 2:column_num]
    end
    fig_width = 100
    if length(ylabel) > 16
        fig_width = 80
    end
    for i in 2:column_num
        plt = lineplot([x for x in data[1]], [y for y in data[i]],
                       title = "$(split(filename,'.')[1])_fig_$(i-1)", name = legend_list[i-1], 
                       xlabel = xlabel, ylabel = ylabel,
                       width = fig_width, height = 20)
        println(plt)
    end
end


main()
