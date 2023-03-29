//Taken from: https://github.com/bokeh/bokeh/blob/branch-3.0/examples/app/export_csv/main.py
//const actual_jitter = jitter_name.properties.value.spec.value
const remove_facet1 = facet1.properties.value.spec.value
const remove_facet2 = facet2.properties.value.spec.value

function table_to_csv(source) {
    const columns = Object.keys(source.data)

    //columns[columns.indexOf("jitter")] = `${actual_jitter}_jitter_points`
    columns[columns.indexOf(`${remove_facet1}_x`)] = remove_facet1
    columns[columns.indexOf(`${remove_facet2}_x`)] = remove_facet2
    columns.splice(columns.indexOf(`${remove_facet1}_y`), 1)
    columns.splice(columns.indexOf(`${remove_facet2}_y`), 1)

    const nrows = source.get_length()
    
    const lines = [columns.join(',')]

    for (let i = 0; i < nrows; i++) {
        let row = [];

        for (let j = 0; j < columns.length; j++) {
            const column = columns[j]

            if(column == remove_facet1) {
                row.push(source.data[`${remove_facet1}_x`][i].toString())

            } else if (column == remove_facet2) {
                row.push(source.data[`${remove_facet2}_x`][i].toString())

            // } else if (column == `${actual_jitter}_jitter_points`) {
            //     row.push(source.data[`jitter`][i].toString())

            } else if (column ==`cat_value`){
                row.push(source.data[column][i].toString().replace(",","|"))

            } else {
                row.push(source.data[column][i].toString())

            }

        }
        lines.push(row.join(','))
    }
    return lines.join('\n').concat('\n')
}

const filename = file_name.properties.text.spec.value
const filetext = table_to_csv(source)
const blob = new Blob([filetext], { type: 'text/csv;charset=utf-8;' })

//addresses IE
if (navigator.msSaveBlob) {
    navigator.msSaveBlob(blob, filename)
} else {
    const link = document.createElement('a')
    link.href = URL.createObjectURL(blob)
    link.download = filename
    link.target = '_blank'
    link.style.visibility = 'hidden'
    link.dispatchEvent(new MouseEvent('click'))
}