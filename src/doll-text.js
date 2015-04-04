var DEBUG=false;

var DRAW_SHIMS=false;

/*
 *
 * Piece-related algorithms and functions.
 *
 * All shim-related dimensions are base-relative (scale: 1 = 1/4").
 *
 */

/** Unit conversions. */
var unitPt = {mm: 72/25.4, cm: 72/2.54, in: 72, pt: 1}; // Conversion from point to unit. 1in = 72pt = 25.4mm

/** SVG attributes for default & bounding box. */
var svgAttr     = {fill: 'none', stroke: 'black', strokeWidth: 0.1},
    svgBBoxAttr = {stroke: 'none'};

/** PDF draw bounding box flag. */
var pdfDrawBBox = false;

/** Maximum seed value. */
var MAX_SEED = 999999;

/** Ratio between shim side (16") and base (1/4"). */
var side = 64;

/** Ratio between shim tip (1/16") and base (1/4"). */
var tip = 0.25;

/** Distance between tip and vanishing point for trapezoidal shim, used to
    compute rotation angle of shims.
    
    tipSide/tip = (tipSide + side)/base, base = 1
    tipSide = tip * (tipSide + side) = tip*tipSide + tip*side
    tipSide - tip*tipSide = tipSide * (1-tip) = tip*side
    tipSide = (tip*side) / (1-tip)
*/
var tipSide = (tip*side) / (1-tip);

/** Angle of shim tips (in radians). Chord is 2*sin(angle/2). */
var shimAngle3 = 2*Math.asin(0.5/side);             /* triangular. */
var shimAngle4 = 2*Math.asin(0.5/(side+tipSide));   /* trapezoidal. */

/** Head-shoulder distance = 5". */
var headShoulder = 20;

/** Shoulder-hip distance = 13 3/8". */
var shoulderHip = 53.5;

/** Hip-calve distance = 14 1/16". */
var hipCalve = 56.25;

/** Bounding box inner margin. */
var bboxMargin = 4;

/** Ten primes used to seed the linear congruential generator. */
var primes = [53, 59, 61, 67, 71, 73, 79, 83, 89, 97];

/** Shuffle factor, improves distribution of LCG on lower digits. Here we used a 
    large prime but any odd number should work given the nature of 2's
    complement arithmetic. */
var shuffle = 57885161;

/** Fonts. Content will be loaded by loadFonts() at startup. */
var fonts = {};

/** Number -> string mapping. Ensure that the array size matches the maximum x value (15). */
var numbers = [
	"ZERO", /* For simplicity */
	"ONE", "TWO", "THREE", "FOUR", "FIVE",
	"SIX", "SEVEN", "EIGHT", "NINE", "TEN",
	"ELEVEN", "TWELVE", "THIRTEEN", "FOURTEEN", "FIFTEEN",
];

/** Number string paths for each font. */
var numberPaths = {};

/**
 *  Load given fonts into global *fonts* array, and compute metrics into *numberMetrics*.
 *  
 *  @param fontList		Font name => url map.
 */
function loadFonts(fontList) {
	$.each(fontList, function(i, url) {
		new FontFile(url, function(font) {
			console.log("Font loaded", url); 
			var name = font.opentype.familyName;
			
			// Sometimes the descender is positive whereas it should be below the baseline at y=0, fix it.
			if (font.opentype.descender > 0) font.opentype.descender = -font.opentype.descender;
			
			fonts[name] = font;
			computeNumberMetrics(name, font);
			
			// Add to select.
			$("#font").append(new Option(name));
		});
	});
}

/**
 *  Compute metrics for all number strings.
 *  
 *  @param name		Font name.
 *  @param font		FontFile object.
 */
function computeNumberMetrics(name, font) {
	numberPaths[name] = [];
	for (var iNumber = 0; iNumber < numbers.length; iNumber++) {
		var glyphs = font.opentype.stringToGlyphs(numbers[iNumber]);
		
		// - Compute raw coords for each glyph.
		var paths = [];
		var maxWidth = 0;
		var lineHeight = font.opentype.ascender-font.opentype.descender;
		var baseline = 0;
		for (var iGlyph = 0; iGlyph < glyphs.length; iGlyph++) {
			var glyph = glyphs[iGlyph];
			
			baseline -= lineHeight; /* Top-down. */
			var commands = [];
			for (var iCommand = 0; iCommand < glyph.path.commands.length; iCommand++) {
				var command = glyph.path.commands[iCommand];
				var coords = [];
				switch (command.type) {
					case 'M': coords = [{x: command.x,  y: (command.y  + baseline)}]; break;
					case 'L': coords = [{x: command.x,  y: (command.y  + baseline)}]; break;
					case 'C': coords = [{x: command.x1, y: (command.y1 + baseline)}, {x: command.x2, y: (command.y2 + baseline)}, {x: command.x, y: (command.y + baseline)}]; break;
					case 'Q': coords = [{x: command.x1, y: (command.y1 + baseline)}, {x: command.x,  y: (command.y  + baseline)}]; break;
				}
				commands[iCommand] = {type: command.type, coords: coords};
			}

			// - Compute real glyph bounding box, overriding values from opentype.js.
			var xMin = Number.POSITIVE_INFINITY, xMax = Number.NEGATIVE_INFINITY,
				yMin = Number.POSITIVE_INFINITY, yMax = Number.NEGATIVE_INFINITY;
			for (var iCommand = 0; iCommand < commands.length; iCommand++) {
				var command = commands[iCommand];
				for (var iCoord = 0; iCoord < command.coords.length; iCoord++) {
					var p = command.coords[iCoord];
					xMin = Math.min(xMin, p.x);
					xMax = Math.max(xMax, p.x);
					yMin = Math.min(yMin, p.y);
					yMax = Math.max(yMax, p.y);
				}
			}
			glyph.xMin = xMin;
			glyph.xMax = xMax;
			glyph.yMin = yMin;
			glyph.yMax = yMax;
			
			var width = glyph.xMax-glyph.xMin;
			maxWidth = Math.max(maxWidth, width);
			
			paths[iGlyph] = commands;
		}
		
		// - Center each glyph horizontally and compute bounding box.
		var xMin = Number.POSITIVE_INFINITY, xMax = Number.NEGATIVE_INFINITY,
			yMin = Number.POSITIVE_INFINITY, yMax = Number.NEGATIVE_INFINITY;
		for (var iGlyph = 0; iGlyph < glyphs.length; iGlyph++) {
			var glyph = glyphs[iGlyph];
			var width = glyph.xMax-glyph.xMin;
			
			var commands = paths[iGlyph];
			for (var iCommand = 0; iCommand < commands.length; iCommand++) {
				var command = commands[iCommand];
				for (var iCoord = 0; iCoord < command.coords.length; iCoord++) {
					var p = command.coords[iCoord];
					
					// - Center each glyph.
					p.x += (maxWidth-width)/2 - glyph.xMin;
					
					// - Compute bounding box.
					xMin = Math.min(xMin, p.x);
					xMax = Math.max(xMax, p.x);
					yMin = Math.min(yMin, p.y);
					yMax = Math.max(yMax, p.y);
				}
			}
		}
		
		// - Normalize coords in (0,0)-(1,1) box.
		for (var iGlyph = 0; iGlyph < glyphs.length; iGlyph++) {
			var glyph = glyphs[iGlyph];
			var width = glyph.xMax-glyph.xMin;
			
			var path = paths[iGlyph];
			for (var iCommand = 0; iCommand < path.length; iCommand++) {
				var command = path[iCommand];
				for (var iCoord = 0; iCoord < command.coords.length; iCoord++) {
					var p = command.coords[iCoord];
					
					p.x = (p.x - xMin) / (xMax - xMin);
					p.y = (p.y - yMin) / (yMax - yMin);
				}
			}
		}
		
		numberPaths[name][iNumber] = paths;
	}
}

/**
 * Generate increment value for linear congruential generator.
 *
 *  @param seed     Seed value whose decimal digits designate prime factors.
 *
 *  @return increment value for lcg()
 */
function lcg_increment(seed) {
    seed %= (MAX_SEED+1);
    if (seed==0) return primes[0];
    
    var c = 1;
    while (seed > 0) {
        c *= primes[seed % 10];
        seed = Math.floor(seed / 10);
    }
    return c;
}

/**
 * Linear congruential generator x_n+1 = (a.x_n + c) mod m.
 *
 * Used to generate a non-repeating sequence of m=x^y integers starting at 0.
 *
 *  @param v    Previous value.
 *  @param c    Increment.
 *  @param x    Number of shims per shim unit.
 *  @param y    Number of shim units/slots per piece.
 *
 *  @return serial number.
 */
function lcg(v, c, x, y) {
    // Number of desired permutations.
    var m = Math.pow(x, y);
    
    // LCG will have a full period if and only if:
    // 1. c and m are relatively prime
    // 2. a-1 is divisible by all prime factors of m
    // 3. a-1 is a multiple of 4 if m is a multiple of 4
    //
    // As m=x^y, prime factors of m are x, if x is prime, or x's prime factors 
    // otherwise.
    // m is multiple of 4 if and only if x is multiple of 2 (assuming y >= 2).
    // #1 is met if x is less than the lowest prime factor used in 
    // lcg_increment().
    
    var a = 2*x+1;  // This guarantees #2 and #3.
    return (a*v+c) % m;
}

/**
 * Generate a random shim permutation.
 *
 *  @param index        Index of piece to generate.
 *  @param c            LCG increment value.
 *  @param x            Number of shims per shim unit.
 *  @param symmetric    Symmetry flag.
 *
 *  @return serial number.
 */
function generatePermutation(index, c, x, symmetric) {
    var y = (symmetric ? 4 : 7);

    // Generate pseudorandom value in [0, 2*max) by calling LCG with sequence
    // number using the previously computed increment value.
    var r = lcg(index * shuffle, c, x, y);
    
    // Digits.
    var digits = "";
    for (var i = 0; i < y; i++) {
        var d = String.fromCharCode(65 + (r % x));
        digits += d;
        if (symmetric && i > 0) digits += d;
        r = Math.floor(r/x);
    }

    return digits;
}

/**
 * Test function validating the LCG-based permutation generator.
 *
 *  Avoid calling with too large values!!!
 *
 *  @param x        Number of shims per shim unit.
 *  @param y        Number of shim units/slots per piece.
 *  @param seed     Seed used to generate LCG increment value.
 */
function testUnicity(x, y, seed) {
    var c = lcg_increment(seed);
    var max = Math.pow(x,y);
    var values = Array();
    var dup = Array();
    for (var i = 0; i < max; i++) {
        var key = lcg(i * shuffle, c, x, y);
        if (typeof(values[key]) === 'undefined') {
            values[key] = [i];
        } else {
            values[key].push(i);
            dup[key] = values[key];
        }
    }
    return {dup: dup, total: Object.keys(values).length};
}

/**
 * Rotate point *p* around center *c* by given *angle*.
 *
 *  @param c        Center of rotation.
 *  @param p        Point to rotate.
 *  @param angle    Rotation angle in radians.
 *
 *  @return Rotated point.
 */
function rotate(c, p, angle) {
    return {
        x: Math.cos(angle) * (p.x-c.x) - Math.sin(angle) * (p.y-c.y) + c.x,
        y: Math.sin(angle) * (p.x-c.x) + Math.cos(angle) * (p.y-c.y) + c.y
    };
}

/**
 * 	Interpolate between *p1* and *p2* at relative position *d*.
 *
 *  @param p1   First point (d=0).
 *  @param p2   Second point (d=1).
 *  @param d    Position.
 *  
 *  @return interpolated point.
 */
function interpolate(p1, p2, d) {
    return {
        x: p1.x + (p2.x-p1.x)*d,
        y: p1.y + (p2.y-p1.y)*d
    };
}

/**
 *  Map point *p* in unit coordinates into trapeze defined by points *p1*-*p2*-*p3*-*p4*.
 *  
 *  @param p1	Top-left (0,0).
 *  @param p2	Bottom-left (0,1).
 *  @param p3	Bottom-right (1,1).
 *  @param p4	Top-right (1,0).
 *  @param p	Source point.
 *  
 *  @return interpolated point.
 */
function trapezeMap(p1, p2, p3, p4, p) {
	// First compute linear interpolations of p.x along (p1,p4) and (p2,p3).
	var p14 = interpolate(p1, p4, p.x),
		p23 = interpolate(p2, p3, p.x);
		
	// Then compute linear interpolation of p.y along (p14,p23).
	return interpolate(p14, p23, p.y);
}

/**
 * Compute a piece from its serial number.
 *
 *  @param sn           The piece serial number.
 *  @param options      Piece options: trapezoidal, font.
 *
 *  @return The piece object.
 */
function computePiece(sn, options) {
    if (sn.length != 7) return;

    //
    // 1. Iterate over slots and build shim coordinates.
    //
    
    var slots = Array(); // Array of slots.
    var angle = options.trapezoidal ? shimAngle4 : shimAngle3;

    // For triangular pieces we need to adjust the shoulder angles so that the 
    // tips don't cross the head's axe of symmetry, else they would overlap.
    // Here we use the law of sins: for a triangle with sides (a,b,c) and 
    // facing angles (A,B,C), a/sinA = b/sinB = c/sinB.
    // A = angle of left head side.
    // a = *side*
    // b = *headShoulder*
    // B = angle to compute = asin(b/a*sinA)
    var headAngleStep = (sn.charCodeAt(0)-64)/2;    // A
    var sideAngleStep = Math.asin(
          headShoulder / side                       // b/a
        * Math.sin(headAngleStep * angle)           // asin(A)
    ) / angle;
    
    for (var iSlot = 0; iSlot < 7; iSlot++) {
        // Left tip corner of first shim.
        var p0_tip = {x: 0, y: side};
            
        // Left base corner of first shim when angle = 0 (vertical).
        var p1_base = {x: 0, y: 0};
        
        // Rotation center.
        var center;
        if (options.trapezoidal) {
            center = {x: 0, y: side+tipSide};
        } else {
            center = p0_tip;
        }

        var shims = Array();
        slots[iSlot] = {shims: shims};
        
        // Iterate over shims.
        var nbShims = sn.charCodeAt(iSlot)-64; /* A=65 */
        var angleStep; // Rotation steps, each of *angle* radians.
        switch (iSlot) {
            case 1:
                // Left shoulder's right side is aligned on head's left side.
                angleStep = -(sn.charCodeAt(0)-64)/2 - nbShims;
                if (!options.trapezoidal) {
                    angleStep += sideAngleStep;
                }
                break;
                
            case 2:
                // Right shoulder's left side is aligned on head's right side.
                angleStep = (sn.charCodeAt(0)-64)/2;
                if (!options.trapezoidal) {
                    angleStep -= sideAngleStep;
                }
                break;
                
            default:
                // Other limbs are hanging symmetrically.
                angleStep= -nbShims/2;
        }               
            
        for (var iShim = 0; iShim < nbShims; iShim++) {
            var p0 = rotate(center, p0_tip, angleStep * angle);
            var p1 = rotate(center, p1_base, angleStep * angle);
            angleStep++;
            var p2 = rotate(center, p1_base, angleStep * angle);
            if (options.trapezoidal) {
                var p3 = rotate(center, p0_tip, angleStep * angle);
                shims[iShim] = [p0, p1, p2, p3];
            } else {
                shims[iShim] = [p0, p1, p2];
            }
        }
    }
    
    
    //
    // 2. Shift slots to junction points.
    //
    
    for (var iSlot = 1; iSlot < 7; iSlot++) {
        var slot = slots[iSlot];
        var prevSlot, prevShim; // Upper shim unit/shim.
        var shiftDistance;      // Distance to interpolate on upper shim unit.
        var junction1;          // Junction point on upper shim unit.
        var junction2;          // Junction point on shim unit to shift.
        var shift;              // Shift vector.
        
        // 2.1. Get previous slot and distance.
        switch (iSlot) {
            case 1:     // Left shoulder.
            case 2:     // Right shoulder.
                prevSlot = slots[0]; // Head.
                shiftDistance = headShoulder/side;
                break;
                
            case 3:     // Left hip.
            case 4:     // Right hip.
                prevSlot = slots[iSlot-2]; // Left/right shoulder.
                shiftDistance = shoulderHip/side;
                break;
                
            case 5:     // Left calve.
            case 6:     // Right calve.
                prevSlot = slots[iSlot-2]; // Left/right hip.
                shiftDistance = hipCalve/side;
                break;
        }
        
        // - 2.2. Compute junction points.
        switch (iSlot) {
            case 1:     // Left shoulder.
            case 3:     // Left hip.
            case 5:     // Left calve.
                // Junction point on left side of upper shim unit.
                prevShim = prevSlot.shims[0];
                junction1 = interpolate(
                    prevShim[1], 
                    prevShim[0], 
                    shiftDistance
                );
                
                // Junction point = right base corner.
                junction2 = slot.shims[slot.shims.length-1][2];
                break;
                
            case 2:     // Right shoulder.
            case 4:     // Right hip.
            case 6:     // Right calve.
                // Junction point on right side of upper shim unit.
                prevShim = prevSlot.shims[prevSlot.shims.length-1];
                junction1 = interpolate(
                    prevShim[2], 
                    prevShim[3%prevShim.length], 
                    shiftDistance
                );
                
                // Junction point = left base corner.
                junction2 = slot.shims[0][1];
                break;
        }
        
        // 2.3. Connect junction points by shifting slots.
        shift = {
            x: junction1.x - junction2.x,
            y: junction1.y - junction2.y,
        };
        
        for (var iShim = 0; iShim < slot.shims.length; iShim++) {
            var shim = slot.shims[iShim];
            for (i = 0; i < shim.length; i++) {
                shim[i].x += shift.x;
                shim[i].y += shift.y;
            }
        }
    }
    
    
    //
    // 3. Compute bounding box.
    //
    
    var x  = Number.POSITIVE_INFINITY, y  = Number.POSITIVE_INFINITY, 
        x2 = Number.NEGATIVE_INFINITY, y2 = Number.NEGATIVE_INFINITY;
    for (var iSlot = 0; iSlot < slots.length; iSlot++) {
        var slot = slots[iSlot];
        for (var iShim = 0; iShim < slot.shims.length; iShim++) {
            var shim = slot.shims[iShim];
            for (i = 0; i < shim.length; i++) {
                x = Math.min(x, shim[i].x);
                y = Math.min(y, shim[i].y);
                x2 = Math.max(x2, shim[i].x);
                y2 = Math.max(y2, shim[i].y);
            }
        }
    }
    
    // Inner margin.
    x  -= bboxMargin; y  -= bboxMargin;
    x2 += bboxMargin; y2 += bboxMargin;
    
    
    //
    // 4. Ensure that all coords are positive. Other code portions depend on it
    //    (e.g. HTML SVG objects & PDF layout).
    //
    
    for (var iSlot = 0; iSlot < slots.length; iSlot++) {
        var slot = slots[iSlot];
        for (var iShim = 0; iShim < slot.shims.length; iShim++) {
            var shim = slot.shims[iShim];
            for (i = 0; i < shim.length; i++) {
                shim[i].x -= x;
                shim[i].y -= y;
            }
        }
    }
    x2 -= x; x = 0;
    y2 -= y; y = 0;
    
	//
	// 5. Fit text into slots.
	//
	
    for (var iSlot = 0; iSlot < slots.length; iSlot++) {
        var slot = slots[iSlot];
		
		// Output trapeze coordinates (or triangle when both ends are the same).
		var firstShim = slot.shims[0], lastShim = slot.shims[slot.shims.length-1];
		var p1 = firstShim[0],
			p2 = firstShim[1],
			p3 = lastShim[2],
			p4 = lastShim[options.trapezoidal ? 3 : 0];
		
		// Iterate over number string's glyph paths.
		var glyphs = numberPaths[options.font][slot.shims.length];
		slot.glyphs = [];
		for (var iGlyph = 0; iGlyph < glyphs.length; iGlyph++) {
			var glyph = glyphs[iGlyph];
			var segments = [];
			
			// Iterate over path commands.
			for (var iCommand = 0; iCommand < glyph.length; iCommand++) {
				var command = glyph[iCommand];
				
				// Map segment coordinates to output trapeze.
				var segment = {type: command.type, coords: []};
				for (var iCoord = 0; iCoord < command.coords.length; iCoord++) {
					segment.coords[iCoord] = trapezeMap(p1, p2, p3, p4, command.coords[iCoord]);
				}
				segments[iCommand] = segment;
			}
			slot.glyphs[iGlyph] = segments;
		}
	}
        
    //
    // Done!
    //
    
    return {sn: sn, slots: slots, bbox: {x: x, y: y, x2: x2, y2: y2}};
}

/**
 * Output a piece as SVG.
 *
 *  @param piece        The piece data.
 *  @param element      DOM element for output (optional).
 *
 *  @return Snap object.
 */
function drawSVG(piece, element) {
    var svg = Snap(element);
    svg.clear();
    svg.attr(svgAttr);
    for (var iSlot = 0; iSlot < piece.slots.length; iSlot++) {
        var slot = piece.slots[iSlot];
		
		if (DRAW_SHIMS) {
			for (var iShim = 0; iShim < slot.shims.length; iShim++) {
				var shim = slot.shims[iShim];
				var coords = Array();
				for (var i = 0; i < shim.length; i++) {
					coords.push(shim[i].x, shim[i].y);
				}
				svg.polygon(coords);
			}
		}
		
		// Build SVG path for each glyph.
		for (var iGlyph = 0; iGlyph < slot.glyphs.length; iGlyph++) {
			var glyph = slot.glyphs[iGlyph];
			var path = "";
			for (var iSegment = 0; iSegment < glyph.length; iSegment++) {
				var segment = glyph[iSegment];
				path += segment.type;
				for (var iCoord = 0; iCoord < segment.coords.length; iCoord++) {
					var c = segment.coords[iCoord];
					path += " " + c.x + " " + c.y;
				}
			}
			svg.path(path).attr({fill: 'black', stroke: 'none'});
		}
    }
    svg.rect(
        piece.bbox.x, piece.bbox.y,
        piece.bbox.x2-piece.bbox.x,
        piece.bbox.y2-piece.bbox.y
    ).attr(svgBBoxAttr);
    return svg;
}

/**
 * Draw a piece into a PDF document.
 *
 *  @param piece        The piece data.
 *  @param pdf          jsPDF document.
 *  @param scale        Scaling factor.
 *  @param offX, offY   Position of top-left corner.
 */
function drawPDF(piece, pdf, scale, offX, offY) {
    // Line width. Use same for shims and bbox.
    pdf.setLineWidth(0.05*scale);
    
    for (var iSlot = 0; iSlot < piece.slots.length; iSlot++) {
        var slot = piece.slots[iSlot];
		
		if (DRAW_SHIMS) {
			for (var iShim = 0; iShim < slot.shims.length; iShim++) {
				var shim = slot.shims[iShim];
				var lines = Array();
				for (var i = 0; i < shim.length; i++) {
					lines.push([
						shim[(i+1)%shim.length].x-shim[i].x,
						shim[(i+1)%shim.length].y-shim[i].y
					]);
				}
				pdf.lines(
					lines,
					shim[0].x*scale+offX, shim[0].y*scale+offY,
					[scale, scale],
					'D'
				);
				
			}
		}

		// Build PDF path for each glyph.
		for (var iGlyph = 0; iGlyph < slot.glyphs.length; iGlyph++) {
			pdf.path(slot.glyphs[iGlyph], offX, offY, [scale, scale], 'F');
		}
    }
	
    if (pdfDrawBBox) {
        pdf.rect(
            piece.bbox.x*scale+offX, piece.bbox.y*scale+offY, 
            (piece.bbox.x2-piece.bbox.x)*scale, (piece.bbox.y2-piece.bbox.y)*scale, 
            'D'
        );
    }
}

/**
 *  Extends jsPDF with path function.
 *  
 *  @param segments		Path segments.
 *  @param x, y			Offset from origin.
 *  @param scale		[x, y] scaling factor.
 * 	@param style		jsPDF style argument.
 *
 *  @return this
 */
jsPDF.API.path = function(segments, x, y, scale, style) {
	style = this.internal.getStyle(style);
	scale = typeof(scale) === 'undefined' ? [1, 1] : scale;
	
	var scalex = scale[0];
	var scaley = scale[1];
	
	var out = this.internal.write;
	var f3 = function (number) {
		return number.toFixed(3);
	};
	var k = this.internal.scaleFactor;
	var pageHeight = this.internal.pageSize.height;
	
	for (var iSegment = 0; iSegment < segments.length; iSegment++) {
		var segment = segments[iSegment];
		for (var iCoord = 0; iCoord < segment.coords.length; iCoord++) {
			var c = segment.coords[iCoord];
			out(f3((c.x * scalex + x) * k) + ' ' + f3((pageHeight - (c.y * scaley + y)) * k) + ' ');
		}
		switch (segment.type) {
			case 'M': out('m'); break;
			case 'L': out('l'); break;
			case 'C': out('c'); break;
			case 'Q': out('v'); break;
			case 'Z': out('h'); break;
		}
	}
	out(style);
	return this;
}

/**
 * Generate a multi-page PDF from a set of pieces.
 *
 *
 *  @param pieceOptions     Piece options: trapezoidal, font.
 *  @param printOptions     Print options:
 *                          - orient    Orientation ('portrait', 'landscape').
 *                          - format    Page format ('a3', 'a4','a5' ,'letter' ,'legal').
 *                          - sides     Print mode ('single', 'double')
 *                          - margins   Margins in unit values {top, bottom, left, right}
 *                          - padding   Padding between pieces in unit values
 *                          - unit      Base measurement unit ('mm', 'cm', 'in', 'pt')
 *                          - justif    Justification ('left', 'center', 'right')
 *                          - cols      Minimum number of columns per page.
 *                          - rows      Minimum number of rows per page.
 *                          - compoPos  Composition number position ('none', 'header','footer').
 *                          - pageNbPos Page number position ('none', 'header','footer').
 *                          - labelPos  Piece S/N label position ('none', 'top','bottom').
 *  @param limits           Output limits:
 *                          - maxPieces        Maximum overall number of pieces to print.
 *                          - maxPiecesPerDoc  Maximum number of pieces per document.
 *                          - maxPagesPerDoc   Maximum number of pages per document.
 *  @param onprogress       Progress callback, called with args (nb, nbPrint, page, nbPages, doc, nbDocs).
 *  @param onfinish         Finish callback.
 */
function piecesToPDF(pieceOptions, printOptions, limits, onprogress, onfinish) {
    var fontSizePt = 10; /* pt */
    
    // Create jsPDF object.
    var pdf = new jsPDF(printOptions.orient, printOptions.unit, printOptions.format);
    var onePt = 1 / pdf.internal.scaleFactor;
    pdf.setFontSize(fontSizePt);
    var fontSizeUnit = fontSizePt * onePt;
    
    // Outer size of the page, i.e. full size minus margins.
    var outerWidth = 
          pdf.internal.pageSize.width 
        - (printOptions.margins.left+printOptions.margins.right);   // Horizontal margins.
    var outerHeight =
          pdf.internal.pageSize.height
        - (printOptions.margins.top+printOptions.margins.bottom);    // Vertical margins.
        
    // Inner size of the page, i.e. outer size minus header and footer.
    var header = (printOptions.compoPos == 'top' || printOptions.pageNbPos == 'top');
    var footer = (printOptions.compoPos == 'bottom' || printOptions.pageNbPos == 'bottom');
    var innerWidth = outerWidth;
    var innerHeight =
          outerHeight
        - (header ? fontSizeUnit+printOptions.padding : 0)          // Header.
        - (footer ? fontSizeUnit+printOptions.padding : 0);         // Footer.
        
    // Compute scaling and actual number of rows/cols.
    var availWidth = innerWidth - printOptions.padding*(printOptions.cols-1);
    var availHeight = innerHeight - printOptions.padding*(printOptions.rows-1);
    var pieceWidth = availWidth / printOptions.cols;
    var pieceHeight = availHeight / printOptions.rows - (printOptions.labelPos == 'none' ? 0 : fontSizeUnit);
    var scale = Math.min(pieceWidth/maxWidth, pieceHeight/maxHeight);
    pieceWidth = maxWidth*scale;
    pieceHeight = maxHeight*scale + (printOptions.labelPos == 'none' ? 0 : fontSizeUnit);
    printOptions.cols = Math.floor((innerWidth + printOptions.padding) / (pieceWidth + printOptions.padding));
    printOptions.rows = Math.floor((innerHeight + printOptions.padding) / (pieceHeight + printOptions.padding));
    
    // Max number of pieces per page.
    var nbPiecesPerPage = printOptions.cols*printOptions.rows;
    
    // Actual number of pieces.
    var nbPrint = Math.min(nbSelected, limits.maxPieces);
    
    // Actual number of pages per document.
    var nbPagesPerDoc = Math.min(Math.ceil(limits.maxPiecesPerDoc/nbPiecesPerPage), limits.maxPagesPerDoc);
    
    // Actual number of pages overall.
    var nbPages = Math.ceil(nbPrint/nbPiecesPerPage);
    
    // Actual number of docs.
    var nbDocs = Math.ceil(nbPages/nbPagesPerDoc);
    
    // Variables for output.
    var col = 0, row = 0, nb = 0, page = 1, firstPage = 1, doc = 1;
    
    // Function for header/footer output.
    var compo = x+"-"+(symmetric?"S":"A")+"-"+seed;
    var compoWidth = pdf.getStringUnitWidth(compo) * fontSizeUnit;
    var headerFooter = function() {
        if (DEBUG) {
			pdf.rect(
				printOptions.margins.left, printOptions.margins.top,
				outerWidth, outerHeight,
				'D'
			);
			pdf.rect(
				printOptions.margins.left, printOptions.margins.top + (header ? fontSizeUnit + printOptions.padding : 0),
				innerWidth, innerHeight,
				'D'
			);
		}
        
        var compoX, compoY, compoJustif;
        var pageNbX, pageNbY, pageNbJustif;
        var pagNbWidth = pdf.getStringUnitWidth(page.toString()) * fontSizeUnit;
        
        // Horizontal positions.
        if (printOptions.pageNbPos == printOptions.compoPos && printOptions.compoPos != 'none') {   
            // Composition and page number side-by-side.
            switch (printOptions.justif) {
                case 'center':
                    // Page number on right/outside.
                    if (printOptions.sides == 'double' && (page % 2) == 0) {
                        compoJustif = 'right';
                        pageNbJustif = 'left';
                    } else {
                        compoJustif = 'left';
                        pageNbJustif = 'right';
                    }
                    break;
                
                case 'left': 
                    // Composition on left.
                    compoJustif = 'left';
                    pageNbJustif = 'right';
                    break;
                    
                case 'right': 
                    // Composition on right.
                    compoJustif = 'right';
                    pageNbJustif = 'left';
                    break;
            }
        } else {
            compoJustif = pageNbJustif = printOptions.justif;
        }
        switch (compoJustif) {
            case 'center':
                compoX = printOptions.margins.left + (innerWidth - compoWidth)/2;
                break;
                
            case 'left':
                compoX = printOptions.margins.left;
                break;
                
            case 'right':
                compoX = printOptions.margins.left + (innerWidth - compoWidth);
                break;
        }
        switch (pageNbJustif) {
            case 'center':
                pageNbX = printOptions.margins.left + (innerWidth - pagNbWidth)/2;
                break;
                
            case 'left':
                pageNbX = printOptions.margins.left;
                break;
                
            case 'right':
                pageNbX = printOptions.margins.left + (innerWidth - pagNbWidth);
                break;
        }
        
        // Output composition number.
        switch (printOptions.compoPos) {
            case 'top':
                compoY = printOptions.margins.top + fontSizeUnit - onePt;
                pdf.text(compoX, compoY, compo);
                break;
            case 'bottom':
                compoY = pdf.internal.pageSize.height - printOptions.margins.bottom;
                pdf.text(compoX, compoY, compo);
                break;
        }
        
        // Output page number.
        switch (printOptions.pageNbPos) {
            case 'top':
                pageNbY = printOptions.margins.top + fontSizeUnit - onePt;
                pdf.text(pageNbX, pageNbY, page.toString());
                break;
            case 'bottom':
                pageNbY = pdf.internal.pageSize.height - printOptions.margins.bottom;
                pdf.text(pageNbX, pageNbY, page.toString());
                break;
        }
    };
    
    // Function for drawing a piece given its index.
    var draw = function(i) {
        // Next column.
        if (nb > 0 && ++col >= printOptions.cols) {
            // Next row.
            col = 0;
            if (++row >= printOptions.rows) {
                // Next page.
                row = 0;
                if ((page % nbPagesPerDoc) == 0) {
                    // Next doc.
                    save();
                    pdf = new jsPDF(printOptions.orient, printOptions.unit, printOptions.format);
                    pdf.setFontSize(fontSizePt);
                    page++;
                    firstPage = page;
                } else {
                    pdf.addPage();
                    page++;
                }
                
                // Invert justification on even pages in double-sided mode.
                if (printOptions.sides == 'double') {
                    switch (printOptions.justif) {
                        case 'left':  printOptions.justif = 'right'; break; 
                        case 'right': printOptions.justif = 'left';  break; 
                    }
                }
                headerFooter();
            }
        }
        
        // Compute piece.
        var sn = generatePermutation(i, c, x, symmetric);
        var piece = computePiece(sn, pieceOptions);

        var labelWidth = pdf.getStringUnitWidth(sn) * fontSizeUnit;

        // Offset in gridded layout.
        var offX = printOptions.margins.left + (pieceWidth + printOptions.padding) * col;
        var offY = printOptions.margins.top + (header ? fontSizeUnit + printOptions.padding : 0) + (pieceHeight + printOptions.padding) * row;
        
        // Justification.
        var labelX = offX;
        var shiftRight = innerWidth - (pieceWidth * printOptions.cols) - (printOptions.padding * (printOptions.cols - 1));

        if (DEBUG) {
			switch (printOptions.justif) {
				case 'left':   pdf.rect(offX, offY, pieceWidth, pieceHeight, 'D'); break;
				case 'center': pdf.rect(offX+shiftRight/2, offY, pieceWidth, pieceHeight, 'D'); break;
				case 'right':  pdf.rect(offX+shiftRight, offY, pieceWidth, pieceHeight, 'D'); break;
			}
		}
        
        switch (printOptions.labelPos) {
            case 'top':
                pdf.text(labelX, offY + fontSizeUnit - onePt*2, sn);
                offY += fontSizeUnit;
                break;
            case 'bottom':
                pdf.text(labelX, offY + pieceHeight, sn);
                break;
        }
        drawPDF(piece, pdf, scale, offX, offY);
        nb++;
    }
    
    // Function for periodic saving.
    var save = function() {
        // Save current PDF document.
        saveAs(new Blob([pdf.output()], {type: 'application/pdf'}), compo+"."+firstPage+"-"+page+".pdf");
        onprogress(nb, nbPrint, page, nbPages, doc, nbDocs);
        doc++;
    }
    
    
    //
    // Now output all pieces!
    //
    
    // First page header.
    headerFooter();
    if (defaultSelected) {
        // All pieces but toggled ones.
        var i = 0;
        var step = 100;
        var drawBg = function() {
            for (; i < nbPieces; i++) {
                if (pieceToggle[i]) continue;
                draw(i);

                if (nb >= nbPrint) {
                    save();
                    setTimeout(onfinish, 0);
                    return;
                }
                
                if ((nb % step) == 0) {
                    onprogress(nb, nbPrint, page, nbPages, doc, nbDocs);
                    setTimeout(drawBg, 0);
                    i++;
                    break;
                }
            }
        }
        drawBg();
    } else {
        // Only toggled pieces
        // Note: no progress here, we expect the number of pieces to be small.
        for (var i in pieceToggle) {
            draw(parseInt(i));

            if (nb >= nbPrint) {
                save();
                setTimeout(onfinish, 0);
                return;
            }
        }
    }
}

/**
 * Generate a Zip archive of SVG files from a set of pieces.
 *
 *  @param pieceOptions     Piece options: trapezoidal.
 *  @param limits           Output limits:
 *                          - maxPieces        Maximum overall number of pieces to export.
 *                          - maxPiecesPerZip  Maximum number of pieces per Zip file.
 *  @param onprogress       Progress callback, called with args (nb, nbPrint, page, nbPages, doc, nbDocs).
 *  @param onfinish         Finish callback.
 */
function piecesToZip(pieceOptions, limits, onprogress, onfinish) {
    // Create JSZip object.
    var zip = new JSZip();
    
    // Actual number of pieces.
    var nbSvg = Math.min(nbSelected, limits.maxPieces);
    
    // Actual number of Zip files.
    var nbFiles = Math.ceil(nbSvg/limits.maxPiecesPerZip);
    
    // Output each piece as SVG file.
    var compo = x+"-"+(symmetric?"S":"A")+"-"+seed;
    var svgTmp = $("#tmpSvg svg")[0];
    var nb = 0, file = 1;
    var generateSvg = function(i) {
        if (nb > 0 && (nb % limits.maxPiecesPerZip) == 0) {
            // Next file.
            save();
            zip = new JSZip();
        }
        
        // Generate SVG from piece.
        var sn = generatePermutation(i, c, x, symmetric);
        var piece = computePiece(sn, pieceOptions);
        var svg = drawSVG(piece, svgTmp);
        svg.attr('viewBox', 
                    piece.bbox.x 
            + " " + piece.bbox.y 
            + " " + (piece.bbox.x2-piece.bbox.x) 
            + " " + (piece.bbox.y2-piece.bbox.y)
        );
        
        // Add SVG to Zip file.
        zip.file(sn + ".svg", svg.outerSVG());
        nb++;
    }
    var save = function() {
        saveAs(zip.generate({type: 'blob', compression: 'DEFLATE'}), compo+((nbFiles > 1) ? "."+file : "")+".zip");
        onprogress(nb, nbSvg, undefined, undefined, file, nbFiles);
        file++;
    }
    
    if (defaultSelected) {
        // All pieces but toggled ones.
        var i = 0;
        var step = 100;
        var generateBg = function() {
            for (; i < nbPieces; i++) {
                if (pieceToggle[i]) continue;
                generateSvg(i);

                if (nb >= nbSvg) {
                    save();
                    setTimeout(onfinish, 0);
                    return;
                }
                
                if ((nb % step) == 0) {
                    onprogress(nb, nbSvg, undefined, undefined, file, nbFiles);
                    setTimeout(generateBg, 0);
                    i++;
                    break;
                }
            }
        }
        generateBg();
    } else {
        // Only toggled pieces
        // Note: no progress here, we expect the number of pieces to be small.
        for (var i in pieceToggle) {
            generateSvg(parseInt(i));

            if (nb >= nbSvg) {
                save();
                setTimeout(onfinish, 0);
                return;
            }
        }
    }
}


/*
 *
 * Interface functions.
 *
 */

/** Handles. */
var x, y, symmetric;

/** Maximum theoretical piece width/height. */
var maxWidth, maxHeight;

/** Number of generated pieces. */
var nbPieces;

/** Permutation seed. */
var seed;

/** LCG increment value generated from seed. */
var c;

/** Columns and rows to display. */
var columns, rows;

/** Column class for piece elements. */
var colClass;

/** Paging. */
var nbPages, nbPerPage;

/** Default selection state. */
var defaultSelected;

/** Piecewise selection toggle state. */
var pieceToggle;

/** Number of toggled pieces. We can't rely on pieceToggle.length because JS 
 *  may switch between vector and object for sparse array storage. */
var nbToggle;

/** Number of selected pieces. */
var nbSelected;

/**
 * Validation for number inputs. Replace the input value with a reasonable
 * number:
 *
 *  - floating point strings are rounded to the nearest step value toward zero.
 *  - final value is kept between min and max.
 *  - non-numeric strings are replaced by zero or the min value.
 */
function validateNumber() {
    var min = $(this).attr('min'); if (!$.isNumeric(min)) min = Number.NEGATIVE_INFINITY;
    var max = $(this).attr('max'); if (!$.isNumeric(max)) max = Number.POSITIVE_INFINITY;
    var step = $(this).attr('step'); if (!$.isNumeric(step)) step = 1;
    this.value = parseFloat((Math.round(parseFloat(this.value)/step)*step).toPrecision(12));
    this.value = Math.max(min, Math.min(max, this.value));
}

/**
 * Ensure that permutation is not too large. Else disable interface elements.
 */ 
function validatePermutationSize() {
    var x = parseInt($("#x").val());
    var symmetric = $("#symmetric").prop('checked');
    var y = (symmetric?4:7);
    var nbPieces = Math.pow(x,y);
    if (nbPieces > Number.MAX_SAFE_INTEGER) {
        // Permutation too large.
        $("#generate").removeClass("btn-default").addClass("btn-danger").prop('disabled', true);
        $("#x, #symmetric").parent().addClass("has-error bg-danger");
        $("#message").addClass("panel-body").html("<div class='alert alert-danger'><span class='glyphicon glyphicon-warning-sign'></span> Permutation size too large!</div>");
    } else {
        $("#generate").removeClass("btn-danger").addClass("btn-primary").prop('disabled', false);
        $("#x, #symmetric").parent().removeClass("has-error bg-danger");
        $("#message").removeClass("panel-body").empty();
    }
}

/**
 * Generate a new set of pieces.
 */
function generatePieces() {
    $("#zip").prop('disabled', false);
    $("#print").prop('disabled', false);

    // Get algorithm handles.
    x = parseInt($("#x").val());
    symmetric = $("#symmetric").prop('checked');
    y = (symmetric?4:7);

    // Number of pieces to generate.
    nbPieces = parseInt($("#nbPieces").val());
    if ($("#max").prop('checked')) {
        // Use max number of pieces.
        nbPieces = Math.pow(x,y);
    }
    
    // Maximum theoretical piece width/height.
    maxWidth = x*7 + 2*bboxMargin;//FIXME better formula?
    maxHeight = headShoulder+shoulderHip+hipCalve+side + 2*bboxMargin;

    // Get/generate seed.
    if ($("#random").prop('checked')) {
        // Generate random seed.
        seed = Math.floor(Math.random() * (MAX_SEED+1));
        $("#seed").val(seed);
    }
    seed = parseInt($("#seed").val());
    seed %= (MAX_SEED+1);
    
    // LCG increment value.
    c = lcg_increment(seed);
    
    // Set default selection state.
    defaultSelected = true;
    pieceToggle = Array();
    nbToggle = 0;
    updateSelected();
    
    // Adjust column layout.
    columns = parseInt($("#columns").val());
    switch (columns) {
        case 0:
            // Automatic, use 4-column responsive layout.
            columns = 4;
            colClass = "col-xs-12";
            if (nbPieces >= 2) {
                colClass += " col-sm-6";
            }
            if (nbPieces >= 3) {
                colClass += " col-md-4";
            }
            if (nbPieces >= 4) {
                colClass += " col-lg-3";
            }
            break;
            
        case 1:
            colClass = "col-xs-12";
            break;
            
        case 2:
            colClass = "col-xs-12 col-sm-6";
            break;
            
        case 3:
            colClass = "col-xs-12 col-sm-4";
            break;
            
        case 4:
            colClass = "col-xs-12 col-sm-3";
            break;
            
        case 6:
            colClass = "col-xs-12 col-sm-2";
            break;
    }
    rows = parseInt($("#rows").val());
    
    // Paging.
    nbPerPage = columns*rows;
    nbPages = Math.ceil(nbPieces/nbPerPage);

    // Display first page
    displayPieces(0);
}

/** 
 * Display pieces for a given page.
 *
 *  @param page     Page number (zero-indexed).
 */
function displayPieces(page) {
    // Sanity check.
    page = Math.max(0, Math.min(page, nbPages-1));
    
    // Display toolbar.
    $("#toolbar").removeClass("hidden");
    
    // Display pager.
    var $pager = $("#pager");
    $pager.empty();
    if (nbPages > 1) {
        var pager = "<ul class='pagination pagination-sm'>";
        pager += "<li" + (page==0?" class='disabled'":"") + "><a href='javascript:displayPieces(" + Math.max(0,page-1) + ")'>&laquo;</a></li>";
        for (var i = 0; i < nbPages; i++) {
            if (nbPages > 10) {
                // Limit buttons to 10, add ellipses for missing buttons.
                if (page < 5) {
                    if (i == 8) {
                        // Ellipsis at end.
                        pager += "<li class='disabled'><span>...</span></li>";
                        i = nbPages-2;
                        continue;
                    }
                } else if (page >= nbPages-5) {
                    if (i == 1) {
                        // Ellipsis at beginning.
                        pager += "<li class='disabled'><span>...</span></li>";
                        i = nbPages-9;
                        continue;
                    }
                } else {
                    if (i == 1) {
                        // Ellipsis at beginning.
                        pager += "<li class='disabled'><span>...</span></li>";
                        i = page-3;
                        continue;
                    } else if (i == page+3) {
                        // Ellipsis at end.
                        pager += "<li class='disabled'><span>...</span></li>";
                        i = nbPages-2;
                        continue;
                    }
                }
            }
            pager += "<li" + (i==page?" class='active'":"") + "><a href='javascript:displayPieces(" + i + ")'>" + (i+1) + "</a></li>";
        }
        pager += "<li" + (page==nbPages-1?" class='disabled'":"") + "><a href='javascript:displayPieces(" + Math.min(page+1,nbPages-1) + ")'>&raquo;</a></li>";
        pager += "</ul>";
        $pager.append(pager);
    }    
    
    // Clear existing pieces.
    var $pieces = $("#pieces");
    $pieces.empty();
    
    // Generate piece output elements.
    var begin = nbPerPage*page;
    var end = Math.min(begin+nbPerPage, nbPieces);
    for (var i = begin; i < end; i++) {
        // Selection state.
        var selected = defaultSelected;
        if (pieceToggle[i]) selected = !selected;
        
        var piece = "<div id='piece-" + i + "' class='form-inline piece " + colClass + "'>";
        piece += "<div class='thumbnail'>";
        piece += "<input id='piece-select-" + i + "' class='piece-select' data-piece='" + i + "' type='checkbox' onclick='togglePiece(" + i + ")' " + (selected?" checked":"") + "/> ";
        piece += "<label for='piece-select-" + i + "'>";
        piece += "<svg xmlns='http://www.w3.org/2000/svg' version='1.1'></svg><br/>";
        piece += "</label>";
        piece += "<div class='input-group input-group-sm'>";
        piece += "<input type='text' class='form-control sn' readonly placeholder='Piece S/N' value='" + generatePermutation(i, c, x, symmetric) + "' size='7'/>";
        piece += "<span class='input-group-btn'><button type='button' class='btn btn-default' onclick='downloadSVG($(this).parent().parent().find(\".sn\").val().trim())'><span class='glyphicon glyphicon-download'></span> SVG</button></span>"
        piece += "</div>";
        piece += "</div>";
        piece += "</div>";
        $pieces.append(piece);
    }
    
    // Display pieces.
    $pieces.find(".piece").each(function(index, element) {
        updatePiece(element);
    });
}

/**
 * Update existing piece when some parameter changes.
 */
function updatePieces() {
    $("#pieces .piece").each(function(index, element) {
        updatePiece(element);
    });
}

/**
 * Compute & output piece from its S/N.
 *
 *  @param element    The containing element.
 */
function updatePiece(element) {
    var sn = $(element).find(".sn").val().trim();
    
    // Generate piece.
    var piece = computePiece(sn, {
        trapezoidal: $("#trapezoidal").prop('checked'),
		font: 		 $("#font").val(),		
    });
    
    // Output to SVG.
    var svg = drawSVG(piece, $(element).find("svg")[0]);
    
    // Adjust viewbox so that all pieces are centered and use the same scale.
    svg.attr('viewBox', 
                ((piece.bbox.x2-piece.bbox.x)-maxWidth)/2
        + " " + ((piece.bbox.y2-piece.bbox.y)-maxHeight)/2
        + " " + maxWidth 
        + " " + maxHeight);
}

/**
 * Toggle select state of given piece.
 *
 *  @param piece    Piece number to toggle.
 */
function togglePiece(piece) {
    if (pieceToggle[piece]) {
        delete pieceToggle[piece];
        nbToggle--;
    } else {
        pieceToggle[piece] = true;
        nbToggle++;
    }
    updateSelected();
}

/**
 * Check/uncheck visible pieces.
 *
 *  @param  check   Whether to check or uncheck pieces.
 */
function checkVisible(check) {
    $(".piece-select").each(function(index, element) {
        if (element.checked^check) {
            $(element).click();
        }
    });
}

/**
 * Check/uncheck all pieces.
 *
 *  @param  check   Whether to check or uncheck pieces.
 */
function checkAll(check) {
    checkVisible(check);
    defaultSelected = check;
    nbToggle = 0;
    pieceToggle = Array();
    updateSelected();
}

/**
 * Update selected piece counters.
 */
function updateSelected() {
    nbSelected = (defaultSelected ? nbPieces - nbToggle : nbToggle);
    $("#totalPieces").html(nbPieces);
    $("#selectedPieces").html(nbSelected);
    $("#zip").prop('disabled', (nbSelected == 0));
    $("#print").prop('disabled', (nbSelected == 0));
}

/**
 * Download piece as SVG.
 *
 *  @param sn   The piece serial number.
 */
function downloadSVG(sn) {
    // Generate piece.
    var piece = computePiece(sn, {
        trapezoidal: $("#trapezoidal").prop('checked'),
		font: 		 $("#font").val(),		
    });
    
    // Output to SVG.
    var svg = drawSVG(piece, $("#tmpSvg svg")[0]);
    svg.attr('viewBox', 
                piece.bbox.x 
        + " " + piece.bbox.y 
        + " " + (piece.bbox.x2-piece.bbox.x) 
        + " " + (piece.bbox.y2-piece.bbox.y)
    );

    blob = new Blob([svg.outerSVG()], {type: "image/svg+xml"});
    saveAs(blob, sn + ".svg");
} 

/**
 * Update progress information during PDF output.
 *
 *  @param ratio    Progress ratio [0,1].
 *  @param piece    Number of pieces output so far.
 *  @param nbPieces Total number of pieces (optional).
 *  @param page     Number of pages output so far (optional).
 *  @param nbPages  Total number of pages (optional).
 *  @param doc      Number of documents output so far (optional).
 *  @param nbDocs   Total number of documents (optional).
 */
function progress(piece, nbPieces, page, nbPages, doc, nbDocs) {
    var percent = (piece/nbPieces)*100;
    $("#progress .progress-bar").attr('aria-valuenow', percent).attr('style','width:'+percent.toFixed(2)+'%').find("span").html(percent.toFixed(0) + "%)");
    $("#progressPiece").html("Piece " + piece + "/" + nbPieces);
    $("#progressPage").html((page && nbPages) ? "Page " + page + "/" + nbPages : "");
    $("#progressDoc").html((doc && nbDocs) ? "Document " + doc + "/" + nbDocs : "");
}

/**
 * Output pieces to PDF. Always use pt units.
 */
function downloadPDF() {
    var units = $("#unit").val();
    
    $("#printDialog").modal('hide');
    $("#progressDialog").modal('show');
    piecesToPDF(
        {
            trapezoidal: $("#trapezoidal").prop('checked'),
			font:		 $("#font").val(),
        },
        {
            orient: $("[name='orient']:checked").val(), 
            format: $("[name='format']:checked").val(),
            sides: $("[name='sides']:checked").val(),

            // Always use pt units, so do conversions upfront.
            unit: 'pt',
            margins: {
                top:    Math.round(parseFloat($("#marginTop").val())    * unitPt[units]),
                bottom: Math.round(parseFloat($("#marginBottom").val()) * unitPt[units]),
                left:   Math.round(parseFloat($("#marginLeft").val())   * unitPt[units]),
                right:  Math.round(parseFloat($("#marginRight").val())  * unitPt[units]),
            },
            padding: Math.round(parseFloat($("#padding").val()) * unitPt[units]),
            
            justif: $("[name='justif']:checked").val(),
            cols: $("#printColumns").val(),
            rows: $("#printRows").val(),
            
            compoPos: $("[name='compoPos']:checked").val(),
            pageNbPos: $("[name='pageNbPos']:checked").val(),
            labelPos: $("[name='labelPos']:checked").val(),
        },
        {
            maxPieces: parseInt($("#maxPieces").val()),
            maxPiecesPerDoc: parseInt($("#maxPiecesPerDoc").val()),
            maxPagesPerDoc: parseInt($("#maxPagesPerDoc").val()),
        },
        progress,
        function() {$("#progressDialog").modal('hide');}
    );
}

/**
 * Output pieces to zipped SVG.
 */
function downloadZip() {
    $("#zipDialog").modal('hide');
    $("#progressDialog").modal('show');
    piecesToZip(
        {
            trapezoidal: $("#trapezoidal").prop('checked')
        },
        {
            maxPieces: parseInt($("#maxZip").val()),
            maxPiecesPerZip: parseInt($("#maxPiecesPerZip").val())
        },
        progress,
        function() {$("#progressDialog").modal('hide');}
    );
}
