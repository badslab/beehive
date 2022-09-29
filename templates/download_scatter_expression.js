//Taken from: https://github.com/bokeh/bokeh/blob/branch-3.0/examples/app/export_csv/main.py
const geneXval = geneX.properties.value.spec.value
const geneYval = geneY.properties.value.spec.value
const geneZval = geneZ.properties.value.spec.value
const numX = num_facet1.properties.value.spec.value
const numY = num_facet2.properties.value.spec.value
const numZ = num_facet3.properties.value.spec.value
const obs_value = obs.properties.value.spec.value
console.log(obs_value,geneXval,geneYval,geneZval,numX,numY,numZ)

function table_to_csv(source) {
    const columns = Object.keys(source.data)
    const nrows = source.get_length()
    columns[columns.indexOf("geneX")] = geneXval
    columns[columns.indexOf("geneY")]= geneYval
    columns[columns.indexOf("geneZ")] = geneZval
    columns[columns.indexOf("num_facetX")] = numX
    columns[columns.indexOf("num_facetY")]= numY
    columns[columns.indexOf("num_facetZ")] = numZ
    columns[columns.indexOf("obs")]= obs_value

    const lines = [columns.join(',')]

    columns[columns.indexOf(geneXval)] = "geneX"
    columns[columns.indexOf(geneYval)]= "geneY"
    columns[columns.indexOf(geneZval)] = "geneZ"
    columns[columns.indexOf(numX)] = "num_facetX"
    columns[columns.indexOf(numY)]= "num_facetY"
    columns[columns.indexOf(numZ)] = "num_facetZ"
    columns[columns.indexOf(obs_value)]= "obs"

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