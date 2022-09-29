//Taken from: https://github.com/bokeh/bokeh/blob/branch-3.0/examples/app/export_csv/main.py

function table_to_csv(source) {
    const columns = Object.keys(source.data)
    const nrows = source.get_length()
    columns.splice(columns.indexOf('padj'), 1)
    columns.splice(columns.indexOf('lfc'), 1)

    const lines = [columns.join(',')]

    for (let i = 0; i < nrows; i++) {
        let row = [];

        for (let j = 0; j < columns.length; j++) {
            const column = columns[j]
                row.push(source.data[column][i].toString())
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