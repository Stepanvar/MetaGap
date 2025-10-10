(function (global) {
    'use strict';

    function onReady(callback) {
        if (document.readyState !== 'loading') {
            callback();
        } else {
            document.addEventListener('DOMContentLoaded', callback);
        }
    }

    function resolveElement(reference) {
        if (!reference) {
            return null;
        }
        if (typeof reference === 'string') {
            return document.querySelector(reference);
        }
        return reference instanceof Element ? reference : null;
    }

    function resolveElements(reference) {
        if (!reference) {
            return [];
        }
        if (typeof reference === 'string') {
            return Array.from(document.querySelectorAll(reference));
        }
        if (reference instanceof Element) {
            return [reference];
        }
        if (typeof reference.length === 'number') {
            return Array.from(reference);
        }
        return [];
    }

    function normaliseText(value) {
        if (value == null) {
            return '';
        }
        return String(value).replace(/\s+/g, ' ').trim();
    }

    function serialiseCell(value, delimiter) {
        var needsQuoting = value.includes('"') || value.includes(delimiter) || value.includes('\n');
        var escaped = value.replace(/"/g, '""');
        if (needsQuoting) {
            return '"' + escaped + '"';
        }
        return escaped;
    }

    function buildFilename(baseName, extension) {
        var name = baseName && baseName.trim() ? baseName.trim() : 'table-export';
        return name + '.' + extension;
    }

    function exportTable(tableApi, format, baseFilename) {
        var delimiter = format === 'tsv' ? '\t' : ',';
        var extension = format === 'tsv' ? 'tsv' : 'csv';
        var headers = tableApi.columns().header().toArray().map(function (cell) {
            return normaliseText(cell && cell.textContent ? cell.textContent : '');
        });
        var rows = tableApi.rows({ search: 'applied' }).nodes().toArray().map(function (row) {
            return Array.from(row.querySelectorAll('td')).map(function (cell) {
                return normaliseText(cell && cell.textContent ? cell.textContent : '');
            });
        });

        var lines = [];
        if (headers.length) {
            lines.push(headers.map(function (value) {
                return serialiseCell(value, delimiter);
            }).join(delimiter));
        }
        rows.forEach(function (row) {
            lines.push(row.map(function (value) {
                return serialiseCell(value, delimiter);
            }).join(delimiter));
        });

        var blob = new Blob([lines.join('\n')], { type: 'text/plain;charset=utf-8' });
        var downloadLink = document.createElement('a');
        downloadLink.href = URL.createObjectURL(blob);
        downloadLink.download = buildFilename(baseFilename, extension);
        document.body.appendChild(downloadLink);
        downloadLink.click();
        document.body.removeChild(downloadLink);
        setTimeout(function () {
            URL.revokeObjectURL(downloadLink.href);
        }, 0);
    }

    function configure(config) {
        if (!global.jQuery || !jQuery.fn || !jQuery.fn.DataTable) {
            return null;
        }

        var tableSelector = config.tableSelector;
        if (!tableSelector) {
            return null;
        }

        var $table = jQuery(tableSelector);
        if (!$table.length) {
            return null;
        }

        var dataTable = $table.hasClass('dataTable')
            ? $table.DataTable()
            : $table.DataTable({
                  dom: config.dom || 'lrtip',
                  pageLength: config.defaultPageLength || 10,
                  order: config.order || [],
              });

        var searchInput = resolveElement(config.searchInput);
        if (searchInput) {
            if (config.searchPlaceholder) {
                searchInput.setAttribute('placeholder', config.searchPlaceholder);
            }
            searchInput.addEventListener('input', function () {
                dataTable.search(this.value).draw();
            });
        }

        var lengthSelect = resolveElement(config.lengthSelect);
        if (lengthSelect) {
            var initialLength = parseInt(lengthSelect.value || dataTable.page.len(), 10);
            if (!Number.isNaN(initialLength)) {
                dataTable.page.len(initialLength).draw();
                lengthSelect.value = String(initialLength);
            }
            lengthSelect.addEventListener('change', function () {
                var value = parseInt(this.value, 10);
                if (!Number.isNaN(value)) {
                    dataTable.page.len(value).draw();
                }
            });
        }

        var exportButtons = resolveElements(config.exportButtons);
        if (exportButtons.length) {
            exportButtons.forEach(function (button) {
                if (button.dataset.tableToolkitBound === 'true') {
                    return;
                }
                button.dataset.tableToolkitBound = 'true';
                button.addEventListener('click', function () {
                    var format = button.getAttribute('data-table-export') || 'csv';
                    var filename = button.getAttribute('data-filename') || config.filename || 'table-export';
                    exportTable(dataTable, format.toLowerCase(), filename);
                });
            });
        }

        return dataTable;
    }

    function init(config) {
        onReady(function () {
            configure(config || {});
        });
    }

    global.TableToolkit = {
        init: init,
    };
})(window);
