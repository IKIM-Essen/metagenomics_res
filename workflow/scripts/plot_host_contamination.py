import pandas as pd
import altair as alt
import sys

## write to log file
sys.stderr = open(snakemake.log[0], "w")

## input file
csv_in = snakemake.input.csv

#input parameter
other_host = snakemake.params.other_host

## output file
host_percentage_html = snakemake.output.html 

## variables
color_red = "#e03e3e"
color_green = "#6aa84f"

## prepare dataframe
df=pd.read_csv(csv_in)

human_cont_df=pd.DataFrame()
human_cont_df["sample"]=df["sample"]
human_cont_df["human"]=df["%human"].divide(100)

if other_host:
    hostname=snakemake.params.hostname
    human_cont_df["host"]=df[f"%{hostname}"].divide(100)


slider = alt.binding_range(
    min=0, max=100, step=0.5, name="maximal accepted host contamination:"
)
# selector = alt.param(name='SelectorName', value=50, bind=slider)
selector = alt.selection_point(
    name="SelectorName", fields=["max_contamination"], bind=slider, value=10
)

# plot for human percentage
if other_host:
    human_base = (
        alt.Chart(human_cont_df)
        .encode(
            alt.X("human:Q")
            .axis(format="%", labelFontSize=12, titleFontSize=12)
            .title("Percentage of human reads")
            .scale(domain=[0, 1]),
            alt.Y("sample:N")
            .axis(labelFontSize=12, titleFontSize=12),
        )
        .add_params(selector)
        .properties(width="container")
        .interactive()
    )

else:
    title="Share of human reads among all reads (after QC)"
    title_object = alt.TitleParams(title, anchor='middle',fontSize=14)

    human_base = (
        alt.Chart(human_cont_df, title=title_object)
        .encode(
            alt.X("human:Q")
            .axis(format="%", labelFontSize=12, titleFontSize=12)
            .title("Percentage of human reads")
            .scale(domain=[0, 1]),
            alt.Y("sample:N")
            .axis(labelFontSize=12, titleFontSize=12),
        )
        .add_params(selector)
        .properties(width="container")
        .interactive()
    )


human_bars = human_base.mark_bar().encode(
    color=alt.condition(
        (alt.datum.human * 100) >= selector.max_contamination,
        alt.value(color_red),
        alt.value(color_green),
    )
)

human_text = human_base.mark_text(
    align="center",
    baseline="middle",
    dx=30,
    fontSize=12,
).encode(
    text=alt.Text("human:Q", format=".2%"),
)

human_full = human_bars + human_text

#if second host
if other_host:
    # plot for other host percentage
    host_base = (
        alt.Chart(human_cont_df)
        .encode(
            alt.X("host:Q")
            .axis(format="%", labelFontSize=12, titleFontSize=12)
            .title(f"Percentage of {hostname} reads")
            .scale(domain=[0, 1]),
            alt.Y("sample:N")
            .axis(labelFontSize=12, titleFontSize=12),
        )
        .add_params(selector)
        .properties(width="container")
        .interactive()
    )

    host_bars = host_base.mark_bar().encode(
        color=alt.condition(
            (alt.datum.host * 100) >= selector.max_contamination,
            alt.value(color_red),
            alt.value(color_green),
        )
    )

    host_text = host_base.mark_text(
        align="center",
        baseline="middle",
        dx=30,
        fontSize=12,
    ).encode(
        text=alt.Text("host:Q", format=".2%"),
    )
    host_full = host_bars + host_text

    title=f"Share of human and {hostname} reads among all reads (after QC)"
    title_object = alt.TitleParams(title, anchor='middle',fontSize=14)

    all_chart=alt.vconcat(host_full, human_full, title=title_object)
    all_chart.save(host_percentage_html)

else:
    human_full.save(host_percentage_html)


