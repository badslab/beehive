
function insertUrlParam(key, value) {
    if (history.pushState) {
        let searchParams = new URLSearchParams(window.location.search);
        searchParams.set(key, value);
        let newurl = window.location.protocol + "//" + window.location.host + window.location.pathname + '?' + searchParams.toString();
        window.history.pushState({path: newurl}, '', newurl);
    }
}

// to remove the specific key
function removeUrlParameter(paramKey) {
    const url = window.location.href
    console.log("url", url)
    var r = new URL(url)
    r.searchParams.delete(paramKey)
    const newUrl = r.href
    console.log("r.href", newUrl)
    window.history.pushState({ path: newUrl }, '', newUrl)
}


function convertToTsv(data, columns) {
    var rv = [ columns.join('\t') ] // header line

    // iterate over the item ids of the first column
    for (i in data[columns[0]]) {
        row = []
        for (k in columns) {
            let val = JSON.stringify(data[columns[k]][i]);
            row.push(val);
        }
        rv.push(row.join('\t'))
    }
    return rv.join("\r\n");
}

function exportToTsv(data, columns, filename)
{
    var tsv = convertToTsv(data, columns);
    var eLink = document.createElement('a');
    eLink.download = filename;
    eLink.style.display = 'none';

    // Converting character content to blob address
    var blob = new Blob([tsv], {type: 'text/csv;charset=utf-8;'});
    eLink.href = URL.createObjectURL(blob);

    // Trigger Click
    document.body.appendChild(eLink);
    eLink.click();

    // remove
    document.body.removeChild(eLink);
};
