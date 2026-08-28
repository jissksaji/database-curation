class WorkflowMain {

    private static final String RESET = '\u001B[0m'

    // xterm-256 neon ramp, hot magenta -> purple -> electric blue -> cyan -> neon green
    private static final List<Integer> NEON = [
        201, 200, 199, 198, 197, 165, 129, 93, 99, 105, 63,
        69, 75, 81, 51, 50, 49, 48, 47, 46, 82, 118,
    ]

    private static final List<String> LOGO = [
        '  _____ _______ _____    ____  _    _ _____    _____  ____    ___  ',
        ' |_   _|__   __/ ____|  / __ \\| |  | |  __ \\  |  __ \\|  _ \\  |__ \\ ',
        '   | |    | | | (___   | |  | | |  | | |__) | | |  | | |_) |    ) |',
        '   | |    | |  \\___ \\  | |  | | |  | |  _  /  | |  | |  _ <    / / ',
        '  _| |_   | |  ____) | | |__| | |__| | | \\ \\  | |__| | |_) |  / /_ ',
        ' |_____|  |_| |_____/   \\____/ \\____/|_|  \\_\\ |_____/|____/  |____|',
    ]

    static void initialise(workflow, params, log) {
        log.info header(workflow, params.monochrome_logs as boolean)
    }

    static String header(workflow, boolean monochrome) {
        """
${logo(monochrome)}
===================================================================
${workflow.manifest.description} | version ${workflow.manifest.version}
===================================================================
"""
    }

    static String logo(boolean monochrome) {
        if (monochrome) {
            return LOGO.join('\n')
        }
        int width = LOGO.collect { it.length() }.max()
        return LOGO.collect { colourise(it, width) }.join('\n')
    }

    // Spread the ramp across the logo width, emitting a code only when the hue changes.
    private static String colourise(String line, int width) {
        StringBuilder out = new StringBuilder()
        int previous = -1
        for (int i = 0; i < line.length(); i++) {
            int colour = NEON[(i * NEON.size()).intdiv(width)]
            if (colour != previous) {
                out << "\u001B[38;5;${colour}m"
                previous = colour
            }
            out << line.charAt(i)
        }
        return out.append(RESET).toString()
    }
}
