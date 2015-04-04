/**
 *  @class FontFile
 *  
 *  Font file class.
 *  
 *  @var url		URL of font file. Must be on same domain.
 *  @var data		Raw binary data (ArrayBuffer).
 *  @var opentype	OpenType info parsed by the OpenType.js library.
 */
 
/**
 * 	FontFile constructor.
 *
 *	@param url			URL of font file. Must be on same domain.
 *	@param onload		Called when data loaded.
 */
function FontFile(url, onload) {
	this.url = url;
	this.data = undefined;
	this.opentype = undefined;
	
	// Load font data.
	this._load(onload);
}

/**
 *  Load font data using AJAX.
 *  
 *	@param onload		Called when data loaded.
 */
FontFile.prototype._load = function(onload) {	
    var request = new XMLHttpRequest();
    request.open('get', this.url, true);
    request.responseType = 'arraybuffer';
	var font = this;
    request.onload = function () {
        if (request.status !== 200) {
            console.log("Font loading failed", font.url, request.statusText);
			return;
		}
		font.data = request.response;
		var info = opentype.parse(font.data);
		if (!info.supported) {
			console.log("Font not supported by OpenType.js", font.url);
		} else {
			font.opentype = info;
		}
		
		// Callback.
		if (onload) onload(font);
    };
    request.send();
	console.log("Font loading scheduled", font.url)
}


/**
 *  Utility callback, schedules a full refresh upon font loading.
 *  
 *  @param font	Loaded FontFile object.
 */
function fontLoaded(font) {
	console.log("Font loaded", font);
	scheduleRefresh();
}