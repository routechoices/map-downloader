const app = require("./api");

const port = 3526;

app.listen(port, () => {
    console.log(`App listening on http://[::]:${port}`)
})